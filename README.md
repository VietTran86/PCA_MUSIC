# Bayesian inference for PCA and MUSIC algorithms with unknown number of sources

Given Y = VA + Z, how to estimate the unknown dimension of V,A optimally, without overfitting ? This is the 50-year-old challenge for popular PCA model.

For the first time, I have found closed-form solution for this challenge via Bayesian method (with linear complexity).

In simulations, we found that SNR = -10 (dB) is the limit of accurate estimation (i.e. non-overfitting) for independent sources. 

In practice, this SNR limit can be estimated from data Y via signal-plus-noise's percentage \tau (i.e. SNR > -10 (dB) <=> \tau < 90%), which means the limit of non-overfitting for independent sources is:

**SNR > -10 (dB)  <=>  "noise's deviation < 3 * source's deviation"**




<p float="left">

  <img src="./2018=PCA MUSIC/Figs/Abstract.png" width="800" />

  <img src="./2018=PCA MUSIC/Figs/Fig1.png" width="400" />
  
  <img src="./2018=PCA MUSIC/Figs/Fig4.png" width="400" />
  
  <img src="./2018=PCA MUSIC/Figs/Fig3.png" width="800" />
  
  <img src="./2018=PCA MUSIC/Figs/Fig2.png" width="400" />
  
  <img src="./2018=PCA MUSIC/Figs/Fig5.png" width="400" />
  
  <img src="./2018=PCA MUSIC/Figs/VietHung_ICASSP_2018.png" width="800" />
  
  <img src="./2018=PCA MUSIC/Figs/Conclusion.png" width="400" />
  
</p>

# Reference:

V.H.Tran and W.Wang, "Bayesian inference for PCA and MUSIC algorithms with unknown number of sources", submitted to IEEE Trans. on Signal Processing, 2018 https://doi.org/10.13140/RG.2.2.21991.50081

V.H.Tran, W.Wang, Y.Luo and J.Chambers, "Bayesian Inference for Multi-Line Spectra in Linear Sensor Array", IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP) 2018, https://ieeexplore.ieee.org/document/8461844
