# Environment
MATLAB R2008a (or higher)
# Add MATLAB Tools
Use ``2DSpecTools`` for **band pass filter**.

Use ``OBNLMpackage`` for **denoising**.

Use ``DIPimage`` for **path opening**.

## Install tools
``2DSpecTools`` and ``OBNLMpackage tools`` have been offered on the GitHub.

Download ``DIPimage`` from http://www.diplib.org/download and install it. 
![alt text](https://github.com/njuRS/picture/blob/master/1549183782(1).jpg?raw=true)

## Set Path
- Click **set path** in Matlab GUI. 
- Add ``2DSpecTools``, ``OBNLMpackage``, ``DIPimage tools``, ``Matlab River Detection code`` folder to **MATLAB environment**. 
- **Save**. 

![alt text](https://github.com/njuRS/picture/blob/master/1549183824(1).jpg?raw=true)

## Open ``Matlab River Detection code`` folder in MATLAB
![alt text](https://github.com/njuRS/picture/blob/master/1549183850(1).jpg?raw=true)

## Change the path of ``DIPimage`` in ``pathopening.m``
![alt text](https://github.com/njuRS/picture/blob/master/1549183861(1).jpg?raw=true)

# Run river detection code
Here, we provide ``run_river_detection.m`` to detect river in single image and ``run_batch_river_detection.m`` for batch detection. 

- Open ``Matlab River Detection code`` folder in MATLAB. 
![alt text](https://github.com/njuRS/picture/blob/master/1549183850(1).jpg?raw=true)
- Open ``run_river_detection.m`` or ``run_batch_river_detection.m``. 
- Write the **image path**. ``test_sentinel2_image.tif`` is provided as test image. 
- Set input **6 parameters**. 
- Click **run**.

If you see these commands in the command window, the river detection code is running successfully.
![alt text](https://github.com/njuRS/picture/blob/master/1549183907(1).jpg?raw=true)
You will get 3 processed images. ``test_sentinel2_image_bandpass_gabor_cpo20.tif`` is the final result.
![alt text](https://github.com/njuRS/picture/blob/master/1549183928(1).jpg?raw=true)
