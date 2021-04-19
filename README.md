# stormsnack
Please do ```git clone https://github.com/QINQINKONG/stormsnack``` to download this repository and it's all good to go.
The following things are contained in this repositary:
- ```cython&dask.ipynb```: The jupyter notebook
- ```integrate.pyx``` : the Cython source file of ```integrate``` function that was used in the notebook
- ```setup_integrate.py```: the ```setup.py``` file for ```integrate.pyx``` source file.
- ```integrate.c```: the C code file cython generated for us.
- ```integrate.cpython-38-x86_64-linux-gnu.so```: extension module that can be imported in python
- ```wetbulb.pyx```: cython source file for wet bulb temperatuure calculation
- ```setup_wetbulb.py```: ```setup.py``` file for wet bulb temperature calculation
- ```wetbulb.c```: the C code file cython generated for uus
- ```wetbulb.cpython-38-x86_64-linux-gnu.so```: extension module that can be imported in python

You can do ```python setup_***.py build_ext --inplace``` to build cython source file which will create a ```.c``` file and ```.so``` file. But since the ```.so``` file has already been included in this repositorary, you don't really need to compile by yourself.
