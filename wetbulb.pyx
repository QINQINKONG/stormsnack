import numpy as np
cimport numpy as np
cimport cython
from libc cimport math
from cython.parallel import prange
from scipy.optimize.cython_optimize cimport brentq
    
cdef double kd,lamda, C, y0, y1, y2
kd = 0.2854
lamda = 3.504
C = 273.15
y0 = 3036.0
y1 = 1.78
y2 = 0.448



cdef double esat(double Tk) nogil:
    return 611.2*math.exp(17.67*(Tk-C)*((Tk-29.65)**(-1)))

cdef double mixrsat(double Tk,double ps) nogil:
    return 0.622*esat(Tk)*((ps - esat(Tk))**(-1))

cdef double vaporpres(double huss, double ps) nogil:
    #huss: specific humidity (kg/kg)
    #ps: surface pressure (Pa)
    #return vapor pressure in Pa
    cdef double r
    r=huss*((1-huss)**(-1))
    return ps*r*((0.622+r)**(-1))

cdef double lcltemp(double Tk,double e) nogil:
    #Tk in K
    #e in Pa
    return 2840.0*(( 3.5*math.log(Tk) - math.log(e/100.0) - 4.805)**(-1)) + 55.0

cdef double thetadl(double Tk, double ps, double e,double Tl,double mixr) nogil:
    return Tk*((100000*((ps-e)**(-1)))**kd)*((Tk*(Tl**(-1)))**(mixr*0.00028))

cdef double thetae(double theta_dl, double Tl, double mixr) nogil:
    return theta_dl*math.exp(((3.036*(Tl**(-1)))-0.00178)*mixr*(1.0 + 0.000448*mixr))

cdef double wb1stguess(double X, double D, double Teq, double ps, double pi) nogil:
    cdef double rs_teq, dlnes_dTeq, wb_temp, k1, k2
    if X > D:
        rs_teq=mixrsat(Teq,ps)
        dlnes_dTeq = 4302.645*((Teq-29.65)**(-2))
        wb_temp = Teq - ((2675.0*rs_teq)*((1.0 + 2675.0*rs_teq*dlnes_dTeq)**(-1)))
    else:
        k1 = pi*(-38.5*pi+137.81)-53.737
        k2 = pi*(-4.392*pi+56.831)-0.384
        if X>=1.0 and X<=D:
            wb_temp = k1-k2*X+C
        elif X>=0.4 and X<1:
            wb_temp = k1-1.21-(k2-1.21)*X+C
        else:
            wb_temp = k1-2.66-(k2-1.21)*X+0.58*(X**(-1))+C
    return wb_temp

cdef double f(double wb, double ps, double rs_wb) nogil:
    cdef double G
    G=(y0*(wb**(-1))-y1)*(rs_wb*(1+y2*rs_wb))
    return ((C*(wb**(-1)))**lamda)*(1.0 - esat(wb)*(ps**(-1)))*math.exp(-lamda*G)

cdef double dfdT(double wb,double ps,double rs_wb) nogil:
    cdef double des_dwb, pminus, rsdT, dGdT
    des_dwb=esat(wb)*4302.645*((wb-29.65)**(-2))
    pminuse = ps - esat(wb) #pminus in Pa
    rsdT=0.622*ps*(pminuse**(-2))*des_dwb
    dGdT = -y0*(rs_wb+y2*rs_wb*rs_wb)*(wb**(-2))+(y0*(wb**(-1))-y1)*(1.0+2.0*y2*rs_wb)*rsdT
    return -lamda*(wb**(-1)+kd*((pminuse)**(-1))*des_dwb+dGdT)*f(wb,ps,rs_wb)
    
ctypedef struct wb_params:
    double C0
    double C1

cdef double fwb(double x, void *args) nogil:
    cdef wb_params *myargs = <wb_params *> args
    cdef double rs,ff,df
    #rs_wb_temp=mixrsat(wb_temp,ps[i,j,k])
    rs=mixrsat(x,myargs.C0)
    ff=f(x,myargs.C0,rs)
    df=dfdT(x,myargs.C0,rs)
    return (ff - myargs.C1)*(df**(-1))

cdef double wb_brentq_wrapper(wb_params args, double xa, double xb, double xtol, double rtol, int mitr) nogil:
    return brentq(fwb, xa, xb, <wb_params *> &args, xtol, rtol, mitr, NULL)


@cython.wraparound(False)
@cython.boundscheck(False)
def wetbulb (double[:,:,::1] Tk,double[:,:,::1] huss, double[:,:,::1] ps, double xtol=0.001, double rtol=0.0, int mitr=1000000):
    cdef Py_ssize_t i, j, k, x_max, y_max, z_max
    x_max = Tk.shape[0]
    y_max = Tk.shape[1]
    z_max = Tk.shape[2]
    result = np.zeros((x_max, y_max, z_max), dtype=np.float64)
    cdef double[:, :, ::1] result_view = result
    cdef double xa,xb,ps_tmp,huss_tmp,Tk_tmp,pi, mixr,e, D,Tl,theta_dl,epott,Teq,X,wb_temp
    cdef wb_params args
    for i in prange(x_max,nogil=True):
        for j in range(y_max):
            for k in range(z_max):
                ps_tmp=ps[i,j,k]
                huss_tmp=huss[i,j,k]
                Tk_tmp=Tk[i,j,k]
                pi = (ps_tmp/100000)**(kd)
                mixr=huss_tmp*((1-huss_tmp)**(-1))*1000 #mixing ratio (g/kg)
                e=vaporpres(huss_tmp,ps_tmp)
                D = (0.1859*ps_tmp/100000 + 0.6512)**(-1)
                Tl = lcltemp(Tk_tmp,e)
                theta_dl=thetadl(Tk_tmp, ps_tmp, e,Tl,mixr)
                epott = thetae(theta_dl,Tl,mixr)
                Teq = epott*pi
                X = (C*(Teq**(-1)))**lamda
                wb_temp=wb1stguess(X, D, Teq,ps_tmp,pi)
                xa=wb_temp-10
                xb=wb_temp+10
                args.C0=ps_tmp
                args.C1=X
                result_view[i,j,k]=wb_brentq_wrapper(args, xa, xb, xtol, rtol, mitr)
    return result
