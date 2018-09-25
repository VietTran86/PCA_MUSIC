%////////////////////////////////////////////

for monte = 1:(setting.MonteCarlo)
    
    Flag.monte.init = (monte == 1);
        
    %////////////////////////////////////////////////////////////////
    
    noise_dB = 10*log10(para.sigma^2); %Noise-to-Signal ratio
    
    data.Y = para.V * para.A + wgn(setting.D,setting.M,noise_dB,'complex');
    
    data.covmat = (data.Y * data.Y');
    
    %////////////////////////////////////////////////////////////////
    
    data.NormY = sum(sum(abs(data.Y).^2));
    
    data.Eigen = sort(eigs(data.covmat,setting.Kmax),'descend');
    
    data.K_AIC = aictest(transpose(data.Y));

 %////////////////////////////////////////////
    
    DTFT = Func_DTFT(setting,para,data);
    
    MUSIC = Func_MUSIC(setting,para,data);
    
end

%//////////////////////////////////////////// result

varNames = fieldnames(DTFT);
for inx_varName = 1:length(varNames)
    varName = varNames{inx_varName};
    
    meanDTFT.(varName)(inxRun) = mean(DTFT.(varName));
    meanMUSIC.(varName)(inxRun) = mean(MUSIC.(varName));
    
end
