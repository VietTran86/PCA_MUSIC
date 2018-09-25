function DTFT = Func_DTFT(in,para,data)
%//////////////////////////////////////////////////////////////////////////

ScanAngles = in.scanAngle:in.scanAngle:180;

wDOA = pi * cos(ScanAngles * pi/180);

%//////////////////////////////////////////////////////////////////////////

V = exp(1i * (1:in.D)' * wDOA);

%------------------------------------

absA = abs(V' * data.Y) / in.D;

PowerA = absA.^2;

SignalPower = sum(PowerA,2)';

%//////////////////////////////////////////////////////////////////////////

[peakPower,inxPeak] = findpeaks(SignalPower);

[peakPowerDescend, inxSpecDescend] = sort(peakPower,'descend');

%//////////////////////////////////////////////////////////////////////////

Kmax = min(length(inxPeak),in.Kmax);

kAngle = ScanAngles(inxPeak(inxSpecDescend(1:Kmax)));

[K_MAP,sigmaKoffset,tauKoffset,logL0_MAP] = func_K_MAP(Kmax,kAngle);

[ERR_Angle_trueK,MSE_Signal_trueK,MSE_SignalTau_trueK] = func_MSE_signal(tauKoffset,kAngle,in.K);
[ERR_Angle_AIC  ,MSE_Signal_AIC  ,MSE_SignalTau_AIC]   = func_MSE_signal(tauKoffset,kAngle,data.K_AIC);
[ERR_Angle_MAP  ,MSE_Signal_MAP  ,MSE_SignalTau_MAP]   = func_MSE_signal(tauKoffset,kAngle,K_MAP);

results();

%//////////////////////////////////////////////////////////////////////////

return

    %//////////////////////////////////////////////////////////////////////

    function results()    

        DTFT.trueK_MSE_Signal    = MSE_Signal_trueK;
        DTFT.trueK_MSE_SignalTau = MSE_SignalTau_trueK;
        DTFT.trueK_ERR_Angle     = ERR_Angle_trueK;
        DTFT.trueK_EST_K         = in.K;
        DTFT.trueK_MSE_sigma     = (sigmaKoffset(in.K+1) - para.sigma)^2;
        DTFT.trueK_tau           = tauKoffset(in.K+1);
        
        DTFT.AIC_MSE_Signal    = MSE_Signal_AIC;
        DTFT.AIC_MSE_SignalTau = MSE_SignalTau_AIC;
        DTFT.AIC_ERR_Angle     = ERR_Angle_AIC;
        DTFT.AIC_EST_K         = data.K_AIC;
        DTFT.AIC_MSE_sigma     = (sigmaKoffset(data.K_AIC+1) - para.sigma)^2;
        DTFT.AIC_tau           = tauKoffset(data.K_AIC+1);
        
        DTFT.MAP_MSE_Signal    = MSE_Signal_MAP;
        DTFT.MAP_MSE_SignalTau = MSE_SignalTau_MAP;
        DTFT.MAP_ERR_Angle     = ERR_Angle_MAP;
        DTFT.MAP_EST_K         = K_MAP;
        DTFT.MAP_MSE_sigma     = (sigmaKoffset(K_MAP+1) - para.sigma)^2;
        DTFT.MAP_tau           = tauKoffset(K_MAP+1);

    end

    %//////////////////////////////////////////////////////////////////////
    
    function [K_MAP,sigmaKoffset,tauKoffset,logL0] = func_K_MAP(Kmax,Angle)
                
        %-----------------------------------------
        
        K0 = 0:Kmax;
        K1 = 1:Kmax;
        
        alphaK1 = K1*in.M;
        
        betaK0 = (in.D - K0)*in.M;
        betaK1 = (in.D - K1)*in.M;
        
        %-----------------------------------------
        
        kNormVA   = zeros(1,Kmax);
        klogProdA = zeros(1,Kmax);
        
        for k=1:Kmax
        
            [V_DOA,A_DOA] = func_VA(Angle(1:k));
            
            kNormVA(k) = sum(sum(abs(V_DOA*A_DOA).^2));
            
            PowerA_DOA = abs(A_DOA).^2;
            
            klogProdA(k) = sum(sum(log(PowerA_DOA)));
        
        end
        
        H0 = data.NormY - [0,kNormVA];
        H1 = H0(2:end);

        logq = log(H1) - log(data.NormY);
        %-----------------------------------------
                
        klogQ1 = zeros(1,Kmax);
        
        for k=1:Kmax
            
            nbeta = 0:(betaK1(k)-1);
                
            nlogNumer = gammaln(nbeta + alphaK1(k));
            
            nlogDenom = gammaln(nbeta+1) + (betaK1(k)-nbeta) * logq(k);
            
            nklogQ = nlogNumer - nlogDenom;
            
            klogQ1(k) = vfunc_logSumExp(nklogQ);
        
        end

        %-----------------------------------------
        
        klogConst = gammaln(betaK1) - gammaln(in.D * in.M);
                
        klogPi = K1 * log(2*pi);
        
        logL1 =  klogConst + klogQ1 - klogPi;

        logL0 = [0,logL1];

        %-----------------------------------------
        
        [~,offsetK_MAP] = max(logL0);
                
        K_MAP = offsetK_MAP-1;
        
        sigmaKoffset = sqrt(H0./(betaK0));
        
        tauKoffset = [1,(H1./betaK1) ./ (kNormVA./alphaK1)];

        %-----------------------------------------
        
    end

    %//////////////////////////////////////////////////////////////////////
    
    function [V,A0] = func_VA(Angle)
        
        wDOA = pi * cos(Angle * pi/180);
        
        V = exp(1i * (1:in.D)' * wDOA);
        
        pseudoV = (V'*V)\V'; 
        
        A0 = pseudoV * data.Y;
        
    end

    %//////////////////////////////////////////////////////////////////////

    function [ErrorRate_Angle,MSE_signal,MSE_signalTau] = func_MSE_signal(tauKoffset,peakAngleDescend,K)
        
        if K == 0
            
            EST_Angles = 0;
            EST_absA = zeros(1,in.M);
            EST_absA_tau = EST_absA;
            
        else
            EST_Angles = peakAngleDescend(1:K);
            [~,EST_A] = func_VA(EST_Angles);
            
            EST_absA = abs(EST_A);
            
            EST_absA_tau = (1-tauKoffset(K+1)) * EST_absA;
        end
        
        MSE_signal = MSEforSignal(EST_Angles,EST_absA,para.sourceAngles,para.A);
        
        MSE_signalTau = MSEforSignal(EST_Angles,EST_absA_tau,para.sourceAngles,para.A);
        
        ErrorRate_Angle = ErrorRateForDOA(EST_Angles,para.sourceAngles,K);
        
    end

    %//////////////////////////////////////////////////////////////////////
    
    function MSE = MSEforSignal(estAngles,estA,trueAngles,trueA)
                   
        [listAngles,inxSort] = sort([estAngles,trueAngles]);
        
        allA = [estA;-trueA];
        allA = allA(inxSort,:);
        
        cumA1 = mean(cumsum(allA).^2,2)';
        
        deltaAngles = [listAngles(2:end),180] - listAngles;
        
        MSE = sum(cumA1 .* deltaAngles) / 180;
        
    end

    %//////////////////////////////////////////////////////////////////////

    function ErrorRate = ErrorRateForDOA(estAngles,trueAngles,K)
        
        if K == 0
            ErrorRate = 1;
        else
            dAngles = abs(trueAngles' - estAngles);
            
            [dAngleMin,inxMin] = min(dAngles);
            
            errClass=zeros(1,in.K);
            
            for k=1:in.K
                
                inxMink = (inxMin == k);
                
                if isempty(find(inxMink,1))
                    errClass(k) = 1;
                else
                    errClass(k) = mean(dAngleMin(inxMink))/180;
                end
            end
            
            ErrorRate = mean(errClass);
        end
    end

    %//////////////////////////////////////////////////////////////////////

    function [LOGX,logMaxX] = vfunc_logSumExp(logX)
        
        logx = logX(:);
        
        logMaxX = max(logx);
        
        if isinf(logMaxX)
            LOGX = logMaxX;
        else
            LOGX = log(sum(exp(logx - logMaxX))) + logMaxX;
        end
        
    end
end