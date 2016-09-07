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
7) figspl_5a.m --> This file generates Figure 5(a) in the mansucript. <br />
8)  figspl_5b.m --> This file generates Figure 5(b) in the mansucript. <br />
9) figspl_5c.m --> This file generates Figure 5(c) in the mansucript. <br />
10)  figspl_5d.m --> This file generates Figure 5(d) in the mansucript. <br />
11)  figspl_5e.m --> This file generates Figure 5(e) in the mansucript. <br />
12)  figspl_5f.m --> This file generates Figure 5(f)  in the mansucript. <br />
13) u_hat_mexhat.m --> This file can generate proposed u_hat approximations for Mexican Hat wavelet with t0=3 and sigma =1 of several orders k1. <br />
14) t_hat_mexhat.m --> This file can generate proposed t_hat approximations for Mexican Hat wavelet with t0=3 and sigma =1 of several orders k1. <br />
15) y_hat_mexhat.m --> This file can generate proposed y_hat approximations for Mexican Hat wavelet with t0=3 and sigma =1 of several orders k1. <br />
16) tou_hat_mexhat.m --> This file can generate proposed tou_hat approximations for Mexican Hat wavelet with t0=3 and sigma =1 of several orders k1. <br />
17) u_hat_gausswav.m --> This file can generate proposed u_hat approximations for Gaussian wavelet with t0=2 and sigma =1 of several orders k1.<br />
18) t_hat_gausswav.m --> This file can generate proposed t_hat approximations for Gaussian wavelet with t0=2 and sigma =1 of several orders k1.<br />
19) y_hat_gausswav.m --> This file can generate proposed y_hat approximations for Gaussian wavelet with t0=2 and sigma =1 of several orders k1.<br />
20) tou_hat_gausswav.m --> This file can generate proposed tou_hat approximations for Gaussian wavelet with t0=2 and sigma =1 of several orders k1.<br />
<br />
<br />
<br />

To generate data points plotted in Figure 6(a) for proposed method in the manuscript run the following files<br />
For order 4 --> Run tou_hat_gausswav.m  after editing k1=4, p=1, beta = -1.02 in the file <br />
For order 5 --> Run tou_hat_gausswav.m after editing k1=5, p=1, beta = -0.96 in the file <br />
For order 6 --> Run y_hat_gausswav.m  after editing k1=6, p=1, beta = -1.09 in the file <br />
For order 7 --> Run y_hat_gausswav.m  after editing k1=7, p=2, beta = -1.88 in the file <br />
For order 8 --> Run y_hat_gausswav.m  after editing k1=8, p=2, beta = -1.9 in the file <br />
The corresponding transfer functions of Pade approach for each order were calculated in MuPad, transferred to MATLAB and then Mean Square Errors were calculated
<br />
<br />

To generate data points plotted in Figure 6(b) for proposed method in the manuscript run the following files<br />
For order 5 --> Run t_hat_mexhat.m  after editing k1=5, p=2, beta = -1.37 in the file <br />
For order 6 --> Run t_hat_mexhat.m  after editing k1=6, p=2, beta = -1.08 in the file <br />
For order 7 --> Run t_hat_mexhat.m  after editing k1=7, p=4, beta = -2.02 in the file <br />
For order 8 --> Run tou_hat_mexhat.m  after editing k1=8, p=4, beta = -3.8 in the file <br />
For order 9 --> Run tou_hat_mexhat.m  after editing k1=9, p=4, beta = -3.83 in the file <br />
The corresponding transfer functions of Pade approach for each order were calculated in MuPad, transferred to MATLAB and then Mean Square Errors were calculated
