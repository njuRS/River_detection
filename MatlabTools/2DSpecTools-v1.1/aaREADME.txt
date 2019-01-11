2DSpecTools for Matlab v1.1
2010/07/06

These Matlab functions perform 2D spectral analysis of gridded topographic data following the procedures described in Perron, J.T., J.W. Kirchner and W.E. Dietrich (2008), Spectral signatures of characteristic spatial scales and non-fractal structure in landscapes, Journal of Geophysical Research - Earth Surface. Please acknowledge the use of this software in any publications by citing this paper.

All code is intended for use with Matlab. The included script mima.m illustrates the use of the included functions. See comments in individual m-files for more details. For help with any function, type 'help [function name]' at the Matlab prompt.

All code copyright (C) 2004-2010 Taylor Perron (perron@mit.edu), except where otherwise indicated. These programs are free software: you can redistribute them and/or modify them under the terms of the GNU General Public License as published by the Free Software Foundation. You should have received a copy of the GNU General Public License along with these programs. If not, see http://www.gnu.org/licenses.


Included files:

Documentation/notes:
  README.txt (this file)

Topographic Data:
  merced.asc       (1)

Matlab Functions:
  bin.m
  CShadePlot.m     (2)
  Detrend.m
  fft2D.m
  fftWin.m
  Hann2D.m
  hillshade.m      (3)
  lsplane.m        (4)
  Make2DFilt.m
  ReadArcGrid.m
  scattercloud.m   (5)
  ShadePlot.m
  SpecFilt2D.m
  SpecPlot1D.m
  SpecPlot2D.m
  WriteArcGrid.m

Matlab Scripts:
  mima.m

Presentations:
  Figures.ppt

Manuscripts:
  Perron-et-al-JGR08.pdf
  Reed-Amundson-2007.pdf


Notes
(1) Data acquired and processed by the National Center for Airborne Laser Mapping (NCALM) at the request of Sarah Reed and Ron Amundson, UC Berkeley. See the included manuscript: Reed, S. and R. Amundson (2007), Sediment, gophers and time:  a model for the origin and persistence of mima mound—vernal pool topography in the Great Central Valley, in R. A. Schlising and D. G. Alexander (Eds), Vernal Pool Landscapes. Studies from the Herbarium, No. 14, California State University, Chico, CA, 15-27.
(2) Based on Mathworks technical note 1215: http://www.mathworks.com/support/tech-notes/1200/1215.shtml
(3) Written by Felix Hebeler. Obtained from the Matlab Central File Exchange:
http://www.mathworks.com/matlabcentral/fileexchange/14863-hillshade
(4) Written by A.B. Forbes and I.M. Smith. Obtained from http://www.eurometros.org/.
(5) Written by Steve Simon. Obtained from the Matlab Central File Exchange:
http://www.mathworks.com/matlabcentral/fileexchange/6037
Distributed under the BSD license (see m-file header for more information).
