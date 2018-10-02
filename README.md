# Bayesian inference for PCA and MUSIC algorithms with unknown number of sources

Given Y = VA + Z, how to estimate the unknown dimension of V, A optimally, without overfitting ? This is the 50-year-old challenge for popular (factor-analysis) PCA model.

For the first time, I have found closed-form solution for this challenge via maximum-a-posterior (MAP) estimate in Bayesian method (i.e. the estimation is fast, with linear complexity). In order to solve this problem, I ended up deriving completely new probability distributions in the Appendix.

In simulations, we found that SNR = -10 (dB) is the limit of accurate estimation (i.e. non-overfitting) for independent sources. 

In practice, this SNR limit can be estimated from data Y via signal-plus-noise's percentage \tau(Y) (i.e. SNR > -10 (dB) <=> \tau(Y) < 90%), which means the limit of non-overfitting for independent sources is:

**SNR > -10 (dB)  <=>  "noise's deviation < 3 * source's deviation"**

P.S: we compared our MAP method with standard MATLAB packages (music and aictest). Everything should be clear in the code. All feedbacks are really welcome!


<p float="left">

  <img src="./2018=PCA MUSIC/Figs/VietHung_PCA_MUSIC.png" width="800" />

  <img src="./2018=PCA MUSIC/Figs/Abstract arxiv.png" width="800" />
  
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

V.H.Tran and W.Wang, "Bayesian inference for PCA and MUSIC algorithms with unknown number of sources", submitted to IEEE Trans. on Signal Processing 2018 https://arxiv.org/abs/1809.10168

V.H.Tran, W.Wang, Y.Luo and J.Chambers, "Bayesian Inference for Multi-Line Spectra in Linear Sensor Array", IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP) 2018, https://ieeexplore.ieee.org/document/8461844

Viet Hung Tran, "Copula Variational Bayes inference via information geometry", IEEE Trans. on information theory 2018 - https://arxiv.org/abs/1803.10998
