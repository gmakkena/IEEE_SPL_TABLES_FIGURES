# IEEE_SPL_TABLES_FIGURES for SPL-20606-2016 	
This document is about the codes that the authors have used to generate results in the current submission titled Nonlinear Sequence Transformation-based Continuous-Time Wavelet Filter Approximation.

MATLAB Version used: MATLAB R2015a <br />
Toolboxes Used by the codes:  <br />
    1) Signal Processing Toolbox  <br />
    2) Control Systems Toolbox  <br />
    3) Symbolic Maths Toolbox <br />
    <br />
    <br />
The following are the MATLAB files that are to be executed to obtain tables and figures in the manuscript.<br />
<br />
1) Table1_col_2_3.m --> This file generates the poles and zeros of the Gaussian wavelet with t0=2, sigma =1 and are displayed in the console. These are found in the coloumns 2 and 3 of the Table 1 in the manuscript. <br />
2) Table1_col_4_5.m --> This file generates the poles and zeros of the Mexican Hat wavelet with t0=2, sigma =1 and are displayed in the console. These are found in the coloumns 4 and 5 of the Table 1 in the manuscript. <br />
3) figspl_1.m --> This file generates the Figure 1 in the manuscript which depicts Absolute percentage change of Taylor coefficents a0 - a5 of Gaussian wavelet as translation is increased. <br />
4) figspl_2.m --> This file generates the Figure 2 in the manuscript which depicts Location of Poles on the complex plane as n is varied for a 5 th order u-transformation of Gaussian wavelet. <br />
5) figspl_3.m --> This file generates Figure 3 in the manuscript which depticts the variation of MSE of the PROPOSED variants as beta is varied. Note that this file takes long time to run because beta value is varied 250 times for each of the four proposed methods. It took about 3min to run this file on a system running on windows 10 with a Intel Xeon processor 3.10GHz and 8Gb RAM<br />
6) figspl_4.m --> This file generates Figure 4 in the manuscript which depticts the variation of Highest Pole quality factor of the PROPOSED variants as beta is varied. Note that this file takes long time to run because beta value is varied 250 times for each of the four proposed methods. It took about 3min to run this file on a system running on windows 10 with a Intel Xeon processor 3.10GHz and 8Gb RAM<br />
7) figspl_5a.m -->
