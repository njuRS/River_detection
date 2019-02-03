
Corresponding Author:

Kang Yang, Xin Lu, Yao Lu

kangyang@nju.edu.cn, xinlu.nju@gmail.com, yaolu.nju@gmail.com

ph: 13814179324

School of Geography and Ocean Science, Nanjing University

Corresponding Author:
Kang Yang, Xin Lu, Yao Lu
kangyang@nju.edu.cn, XXX, XXX
ph: 13814179324
School of Geography and Ocean Science, Nanjing University
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

| parameters | description |
|----|---|
|  |  **WV**: WorldView image  |
|  |  **SPOT**: SPOT image  |
|*sensor*  |  **Sentinel2**: Sentinle-2 image after NDWI calculation  |
|  |  **Landsat**: Landsat panchromatic image  |
|  |  **LandsatNDWILandsat**: image after NDWI calculation  |
|*inverse*  |set 1 to convert dark rivers into bright rivers (e.g., for panchromatic image); set 0 to keep bright rivers (e.g., for NDWI images)|
|*width*  |small river width for Gabor filter, default = 2|
|*ppo_length*  |path opening length, default = 20|
|*histCountThreshold*  |a pixel count threshold to stretch image pixel values, default = 1000|
|*Smooth (optional)*  |used for denoise algorithm. Because denoise algorithm is too slow and scale dependent, is abandoned for large images,  default = 0.7|

**You can add your customized sensor (image) easily. The only requirement is to reset the band pass frequency based on the spatial resolution of your input image.**

- Click **run**.

If you see these commands in the command window, the river detection code is running successfully.
![alt text](https://github.com/njuRS/picture/blob/master/1549183907(1).jpg?raw=true)
You will get 3 processed images. ``test_sentinel2_image_bandpass_gabor_cpo20.tif`` is the final result.
![alt text](https://github.com/njuRS/picture/blob/master/1549183928(1).jpg?raw=true)

For more information, please see the paper: 

Yang, K. , Karlstrom, L., Smith, L.C., Li, M., 2017. Automated high resolution satellite image registration using supraglacial rivers on the Greenland Ice Sheet. IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 10(3): 845-856.

Yang, K. , Li, M., Liu, Y., Cheng, L., Huang, Q., Chen, Y., 2015. River Detection in Remotely Sensed Imagery Using Gabor Filtering and Path Opening. Remote Sensing, 7(7): 8779-8802.

Yang, K. , Li, M., Liu, Y., Cheng, L., Duan, Y., Zhou, M., 2014. River Delineation from Remotely Sensed Imagery Using a Multi-Scale Classification Approach. IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 7(12): 4726-4737.
