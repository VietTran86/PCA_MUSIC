%//////////////////////////////////////////////////////////////////////////
%------------------------------------- % source's direction of angles
para.sourceAngles = setting.offsetAngle + (0:setting.K-1) * setting.interAngle;

para.wDoA = pi * cos(para.sourceAngles * pi/180);
para.V = exp(1i * (1:setting.D)' * para.wDoA);

%------------------------------------- % source's tones

interFFT = setting.offsetFFT + (0:setting.K-1) * setting.interFFT;

bandFFT = interFFT' + (1:setting.bandFFT);

%------------------------------------- % source's amplitudes

para.A = zeros(setting.K,setting.M);

AmplitudeRatio = 1 - (0:setting.K-1) * setting.AmplitudeDecay/setting.K;


for j=1:setting.K
    para.A(j,bandFFT(j,:)) = setting.Amplitude * AmplitudeRatio(j);%*0;
end
%-------------------------------------

para.SNR = 10^(setting.SNR_dB(inxRun)/10);

%------------------------------------- % noise's power

para.SignalPowerDOA = (setting.Amplitude^2 * setting.bandFFT)/setting.M;

para.SignalPower = para.SignalPowerDOA * setting.D;

para.NoisePower = para.SignalPower / para.SNR;

para.sigma  = sqrt(para.NoisePower);

para.sigmaD = para.sigma / sqrt(setting.D);

%-------------------------------------

para.VA = para.V * para.A;
%//////////////////////////////////////////////////////////////////////////


