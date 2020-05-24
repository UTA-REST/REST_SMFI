# REST_SMFI

A set of image analysis functions for identifying differences in bulk and single molecule samples.

### Installing

After you clone the repo the setup.sh script should set up the enviorment.
You will need to run this setup every time you launch a new terminal. 

```bash
source setup.sh
```
If you want to build the Cython functions you will need to go into Core/Cython and source the buils.sh file.
This is needed for the slide histograming function. 


### running

In the Exaples folder is Fluorescence_Trajectories.ipynb, this is an example of analysing an imace stack and is the only full example in the repo the others require your data files. 
