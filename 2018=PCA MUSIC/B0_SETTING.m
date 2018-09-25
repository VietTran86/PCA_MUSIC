%//////////////////////////////////////////////////////////////////////////
%------------------------------------- 
setting.M = 2^12; % No of FFT tones
%------------------------------------- 
setting.D = 100; % No of sensors
%-------------------------------------
setting.K = 5; % true No of sources
%-------------------------------------
setting.Kmax = 10; % maximum No of sources (estimate)
%/////////////////////////////////////////////////////

%/////////////////////////////////////////////////////
setting.Amplitude = 1;
%------------------------------------- 
setting.AmplitudeDecay = 0; % Amplitude's decay rate (\psi = [0,1])
%/////////////////////////////////////////////////////

%/////////////////////////////////////////////////////
setting.scanAngle = 0.1;
%------------------------------------- 
setting.offsetAngle = 10; 
%------------------------------------- % source's angle resolution (degree)
% setting.interAngle = 4;       % correlated DOAs
setting.interAngle = floor((180-setting.offsetAngle)/setting.K);    
%/////////////////////////////////////////////////////

%/////////////////////////////////////////////////////
%------------------------------------- % source's Freq offset (DFT-bin)
setting.offsetFFT = 0;
%-------------------------------------  % No of tones (sources)
setting.bandFFT = floor(setting.M/setting.K);
%------------------------------------- % source's Freq resolution (DFT-bin)
% setting.interFFT = 0;             % complete overlapping (100%)
% setting.interFFT = 1;             % almost overlapping (99.9%)
setting.interFFT = setting.bandFFT; % non-overlapping (0%)
%/////////////////////////////////////////////////////

%/////////////////////////////////////////////////////

setting.MonteCarlo = 1000; % number of Monte Carlo runs

%/////////////////////////////////////////////////////

setting.SNR_dB = -40:5:40; 