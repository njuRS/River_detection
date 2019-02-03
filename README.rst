**Environment**
>MATLAB R2008a (or higher)
Add MATLAB Tools:
Use ``2DSpecTools`` for band pass filter. 
Use OBNLMpackage for denoising. 
Use DIPimage for path opening. 
1.	Install tools: 
2DSpecTools and OBNLMpackage tools have been offered on the GitHub. Download DIPimage from http://www.diplib.org/download and install it. 
 
2.	Set Path: 
1)	Click set path in Matlab GUI. 
2)	Add 2DSpecTools, OBNLMpackage, DIPimage tools, ‘Matlab River Detection code’ folder to MATLAB environment. 
3)	Save. 
 

 
3.	Open ‘Matlab River Detection code’ folder in MATLAB.
 
4.	Change the path of DIPimage in ‘pathopening.m’. 
 
Run river detection code: 
Here, we provide ‘run_river_detection.m’ to detect river in single image and ‘run_batch_river_detection.m’ for batch detection. 
1.	Open ‘Matlab River Detection code’ folder in MATLAB. 
 
2.	Open ‘run_river_detection.m’ or ‘run_batch_river_detection.m’. 
3.	Write the image path. ‘test_sentinel2_image.tif’ is provided as test image. 
4.	Set input 6 parameters. 
1)	sensor: type of your input image, we gave 6 types of input. 
WV: WorldView image; 
SPOT: SPOT image; 
SETSM: ArcticDEM image; 
Sentinel2: Sentinle-2 image after NDWI calculation; 
Landsat: Landsat panchromatic image
LandsatNDWI: Landsat image after NDWI calculation; 
2)	inverse: set 1 to convert dark rivers into bright rivers (e.g., for panchromatic image); set 0 to keep bright rivers (e.g., for NDWI images). 
3)	width: small river width for Gabor filter, default = 2. 
4)	ppo_length: path opening length, default = 20.
5)	histCountThreshold: a pixel count threshold to stretch image pixel values, default = 1000.
6)	Smooth (optional): this parameter is used for denoise algorithm. Because denoise algorithm is too slow and scale dependent, is is abandoned for large images,  default = 0.7
7)	 You can add your customized sensor (image) easily. The only requirement is to reset the band pass frequency based on the spatial resolution of your input image.
5.	Click run. 
If you see these commands in the command window, the river detection code is running successfully. 
 
6.	You will get 3 processed images. ’test_sentinel2_image_bandpass_gabor_cpo20.tif’ is the final result. 
 
