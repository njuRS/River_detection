# Environment
MATLAB R2008a (or higher)
# Add MATLAB Tools
Use ``2DSpecTools`` for **band pass filter**.

Use ``OBNLMpackage`` for **denoising**.

Use ``DIPimage`` for **path opening**.

## Install tools
``2DSpecTools`` and ``OBNLMpackage tools`` have been offered on the GitHub.

Download ``DIPimage`` from http://www.diplib.org/download and install it. 

## Set Path
- Click **set path** in Matlab GUI. 
- Add ``2DSpecTools``, ``OBNLMpackage``, ``DIPimage tools``, ``Matlab River Detection code`` folder to **MATLAB environment**. 
- **Save**. 

## Open ``Matlab River Detection code`` folder in MATLAB

## Change the path of ``DIPimage`` in ``pathopening.m``

# Run river detection code
Here, we provide ``run_river_detection.m`` to detect river in single image and ``run_batch_river_detection.m`` for batch detection. 

- Open ``Matlab River Detection code`` folder in MATLAB. 
- Open ``run_river_detection.m`` or ``run_batch_river_detection.m``. 
- Write the **image path**. ``test_sentinel2_image.tif`` is provided as test image. 
- Set input **6 parameters**. 
- Click **run**.

If you see these commands in the command window, the river detection code is running successfully.

You will get 3 processed images. ``test_sentinel2_image_bandpass_gabor_cpo20.tif`` is the final result. 
