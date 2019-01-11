/* Pierrick Coupe - pierrick.coupe@gmail.com                               */
/* Brain Imaging Center, Montreal Neurological Institute.                  */
/* Mc Gill University                                                      */
/*                                                                         */
/* Copyright (C) 2008 Pierrick Coupe                                       */



/*                 Details on Bayesian NLM filter                         */
/***************************************************************************
 *  The bayesian NLM filter is described in:                               *
 *                                                                         *
 * P. Coupe, P. Hellier, C. Kervrann, C. Barillot.                         *
 * NonLocal Means-based Speckle Filtering for Ultrasound Images.           *
 * IEEE Transactions on Image Processing, 18(10):2221?9, 2009.             *
 ***************************************************************************/


/*                 Details on blockwise NLM filter                        */
/***************************************************************************
 *  The blockwise NLM filter is described in:                              *
 *                                                                         *
 *  P. Coup?, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.     *
 *  An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic*
 *  Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441, * 
 *  Avril 2008                                                             *
 ***************************************************************************/



/*                 This method is patented as follows                     */
/***************************************************************************
* P. Coup?, P. Hellier, C. Kervrann, C. Barillot. Dispositif de traitement *
* d?images am?lior?. INRIA, Patent: 08/02206, 2008.                        *
* Publication No. : WO/2009/133307.                                        *
* International Application No. : PCT/FR2009/000445                        *                                    *
 ***************************************************************************/


Matlab code for Bayesian Non-local means speckle filtering

1) Select the directory where the files have been extracted 
as work directory in matlab.

2) Open the script : SpeckleRemoval.m (adjust the paramters if required)
	
3) Run the script 

4) Use the dialogue interface to select your input image and the ouput directory

5) Check the displayed results

6) The intensity normalized denoised and noisy images should be in the selected output directory


Warning #1: the intensity of input (thus output) image is normalized between 0-1 for writing purpose.
In order to be able to read and write the main image formats and handle the type issues, I decided to normalized the input image to [0-1] for writting. For some formats, Maltab directly write the images between [0-1] to [0-255]. Finally, I deciced to write the normalized input and output images in order to enable comparison and visual assessement of denoising quality. 

Warning #2: this implementation is different to the C++ version used in the TIP paper. 
The original version is licenced and cannot be distributed. 
In consequence, the optimal parameters can differ from the TIP paper.
Tips:
- If you can see structures in the "Residual image" decrease the smoothing parameter h
- If you can see speckle in the "denoised image" increase this parameter


