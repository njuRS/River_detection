% /* Pierrick Coupe - pierrick.coupe@gmail.com                               */
% /* Brain Imaging Center, Montreal Neurological Institute.                  */
% /* Mc Gill University                                                      */
% /*                                                                         */
% /* Copyright (C) 2008 Pierrick Coupe                                       */
%
%
%
% /*                 Details on Bayesian NLM filter                         */
% /***************************************************************************
%  *  The bayesian NLM filter is described in:                               *
%  *                                                                         *
%  * P. Coupe, P. Hellier, C. Kervrann, C. Barillot.                         *
%  * NonLocal Means-based Speckle Filtering for Ultrasound Images.           *
%  * IEEE Transactions on Image Processing, 18(10):2221-9, 2009.             *
%  ***************************************************************************/
%
%
% /*                 Details on blockwise NLM filter                        */
% /***************************************************************************
%  *  The blockwise NLM filter is described in:                              *
%  *                                                                         *
%  *  P. Coupe, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.     *
%  *  An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic*
%  *  Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441, *
%  *  Avril 2008                                                             *
%  ***************************************************************************/
%
%
%
% /*                 This method is patented as follows                     */
% /***************************************************************************
% * P. Coupe, P. Hellier, C. Kervrann, C. Barillot. Dispositif de traitement *
% * d'images ameliore. INRIA, Patent: 08/02206, 2008.                        *
% * Publication No. : WO/2009/133307.                                        *
% * International Application No. : PCT/FR2009/000445                        *
% *                                                                          *
%  ***************************************************************************/

warning off;
clear all;
close all;
clc;

%params
M = 7;      % search area size (2*M + 1)^2
alpha = 3;  % patch size (2*alpha + 1)^2
h = 0.7;    % smoothing parameter [0-infinite].
% If you can see structures in the "Residual image" decrease this parameter
% If you can see speckle in the "denoised image" increase this parameter

offset = 100; % to avoid Nan in Pearson divergence computation
% According to the gain used by your US device this offset can be adjusted.


[namein, pathin] = uigetfile({'*.tif;*.tiff;*.png;*.jpg;*.bmp;*.pgm;*.ppm', 'All Supported Formats'}, 'Open image');
if isequal(namein,0) | isequal(pathin,0)
    disp('User pressed cancel')
else
    
    disp(['Input file : ', fullfile(pathin, namein)])
    
    
    [pathstr, name_s, ext, versn]=fileparts(fullfile(pathin, namein));
    nameout=[name_s '_norm_denoised' ext];
    nameinnorm=[name_s '_norm' ext];
    %namespeckle=[name_s '_speckle' ext];
    
    [nout, pathout] = uiputfile({['*' ext]},'Save as',[pathin nameout]);
    if isequal(nout,0) | isequal(pathout,0)
        disp('User pressed cancel')
    else
        
        disp(['Output file: ', fullfile(pathout, nout)])
        
        img = imread([pathin namein]);
        info = imfinfo([pathin namein]);
        
        [dimxy dimt] = size(size(img));
        if ( dimt > 2)
            img = rgb2gray(img);
        end
        
        % Intensity normalization
        imgd = double(img);
        mini = (min(imgd(:)));
        imgd = (imgd - mini);
        maxi = max(imgd(:));
        imgd = (imgd / maxi) * 255;
        imgd = imgd + offset; % add offset to enable the pearson divergence computation (i.e. avoid division by zero).
        s = size(imgd);
        
        % Padding
        imgd = padarray(imgd,[alpha alpha],'symmetric');
        fimgd=bnlm2D(imgd,M,alpha,h);
        fimgd = fimgd - offset;
        imgd = imgd - offset;
        imgd = imgd(alpha+1: s(1)+alpha, alpha+1: s(2)+alpha);
        fimgd = fimgd(alpha+1: s(1)+alpha, alpha+1: s(2)+alpha);
        
        % Display
        minds = min(imgd(:));
        maxds = max(imgd(:));
        figure;
        imagesc(imgd,[minds maxds]);
        title('Original')
        colormap(gray);
        colorbar;
        figure;
        colormap(gray);
        imagesc(fimgd,[minds maxds]);
        title('Denoised by Bayesian NLM')
        colorbar;
        figure;
        colormap(gray);
        speckle = abs(imgd(:,:) - fimgd(:,:));
        imagesc(speckle);
        title('Residual image')
        colorbar;
        
        
        if (  strcmp(ext , '.gif') || strcmp(ext , '.bmp') || strcmp(ext , '.jpg') || strcmp(ext , '.ppm') ||  strcmp(ext , '.png') || (strcmp(ext , '.pgm') ))
            fimg = uint8(fimgd);
            img = uint8(imgd);
            %speckle = uint8(speckle);
            
        elseif(strcmp(ext , '.tif') || strcmp(ext , '.tiff'))

            fimg = fimgd/ max(imgd(:));
            img = imgd/ max(imgd(:)); 
            %speckle = speckle/ max(imgd(:)); 
        end
        
        imwrite(fimg,[pathout nout]);
        imwrite(img,[pathout nameinnorm]);
        %imwrite(speckle,[pathout namespeckle]);
        
    end
end

