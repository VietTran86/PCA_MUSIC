%////////////////////////////////////////////////////////////////////////////////////////////////
st.Line  = 'LineWidth';
st.Msize = 'MarkerSize';
st.Disp  = 'DisplayName'; 
%--------------------------------------------
spec.DTFT_trueK   = {'g:h',   st.Disp,'DTFT (true K)'     st.Line,1,st.Msize,10};
spec.DTFT_AIC     = {'b--x',   st.Disp,'DTFT (AIC)'    st.Line,1,st.Msize,10};
spec.DTFT_MAP    =  {'k-.+',    st.Disp,'DTFT (MAP)'    st.Line,1,st.Msize,15};
%--------------------------------------------------------------------------
spec.MUSIC_trueK    = {'g:d', st.Disp,'MUSIC (true K)'     st.Line,1,st.Msize,15};
spec.MUSIC_AIC      = {'b--o',    st.Disp,'MUSIC (AIC)'            st.Line,1,st.Msize,13};
spec.MUSIC_MAP      = {'r-.s',   st.Disp,'MUSIC (MAP)'       st.Line,1,st.Msize,17};
%////////////////////////////////////////////////////////////////////////////////////////////////

%//////////////////////////////////////////////////////////////////////////

h.K = figure;
hold on;
%-----------------------------------------
plot(setting.SNR_dB,meanDTFT.trueK_EST_K,spec.DTFT_trueK{:});
plot(setting.SNR_dB,meanDTFT.  AIC_EST_K,spec.DTFT_AIC{:});
plot(setting.SNR_dB,meanDTFT.  MAP_EST_K,spec.DTFT_MAP{:});
plot(setting.SNR_dB,meanMUSIC.trueK_EST_K,spec.MUSIC_trueK{:});
plot(setting.SNR_dB,meanMUSIC.  AIC_EST_K,spec.MUSIC_AIC{:});
plot(setting.SNR_dB,meanMUSIC.  MAP_EST_K,spec.MUSIC_MAP{:});
%-----------------------------------------
xlabel('SNR (dB)');
ylabel('Averaged estimate of number K of sources');

legend
movegui(h.K,'north')

%//////////////////////////////////////////////////////////////////////////

h.DOA = figure;
hold on;
%-----------------------------------------
plot(setting.SNR_dB,meanDTFT.trueK_ERR_Angle ,spec.DTFT_trueK{:});
plot(setting.SNR_dB,meanDTFT.  AIC_ERR_Angle ,spec.DTFT_AIC{:});
plot(setting.SNR_dB,meanDTFT.  MAP_ERR_Angle ,spec.DTFT_MAP{:});
plot(setting.SNR_dB,meanMUSIC.trueK_ERR_Angle ,spec.MUSIC_trueK{:});
plot(setting.SNR_dB,meanMUSIC.  AIC_ERR_Angle ,spec.MUSIC_AIC{:});
plot(setting.SNR_dB,meanMUSIC.  MAP_ERR_Angle ,spec.MUSIC_MAP{:});
%-----------------------------------------
xlabel('SNR (dB)');
ylabel('Error rate of DOA (%)');

legend
movegui(h.DOA,'south')

%//////////////////////////////////////////////////////////////////////////

h.sigma = figure;
hold on;
%-----------------------------------------
semilogy(setting.SNR_dB,sqrt(meanDTFT.trueK_MSE_sigma),spec.DTFT_trueK{:});
semilogy(setting.SNR_dB,sqrt(meanDTFT.  AIC_MSE_sigma),spec.DTFT_AIC{:});
semilogy(setting.SNR_dB,sqrt(meanDTFT.  MAP_MSE_sigma),spec.DTFT_MAP{:});
semilogy(setting.SNR_dB,sqrt(meanMUSIC.trueK_MSE_sigma),spec.MUSIC_trueK{:});
semilogy(setting.SNR_dB,sqrt(meanMUSIC.  AIC_MSE_sigma),spec.MUSIC_AIC{:});
semilogy(setting.SNR_dB,sqrt(meanMUSIC.  MAP_MSE_sigma),spec.MUSIC_MAP{:});
%-----------------------------------------
xlabel('SNR (dB)');
ylabel('RMSE of noise''s deviation \sigma');

set(gca, 'YScale', 'log')

legend
movegui(h.sigma,'northeast')

%//////////////////////////////////////////////////////////////////////////

h.tau = figure;
hold on;
%-----------------------------------------
plot(setting.SNR_dB,meanDTFT.trueK_tau ,spec.DTFT_trueK{:});
plot(setting.SNR_dB,meanDTFT.  AIC_tau ,spec.DTFT_AIC{:});
plot(setting.SNR_dB,meanDTFT.  MAP_tau ,spec.DTFT_MAP{:});
plot(setting.SNR_dB,meanMUSIC.trueK_tau ,spec.MUSIC_trueK{:});
plot(setting.SNR_dB,meanMUSIC.  AIC_tau ,spec.MUSIC_AIC{:});
plot(setting.SNR_dB,meanMUSIC.  MAP_tau ,spec.MUSIC_MAP{:});
%-----------------------------------------
xlabel('SNR (dB)');
ylabel('Noise-to-signal percentage \tau (%)');

legend
movegui(h.tau,'southeast')

%//////////////////////////////////////////////////////////////////////////

h.ML = figure;
hold on;
%-----------------------------------------
semilogy(setting.SNR_dB,sqrt(meanDTFT.trueK_MSE_Signal),spec.DTFT_trueK{:});
semilogy(setting.SNR_dB,sqrt(meanDTFT.  AIC_MSE_Signal),spec.DTFT_AIC{:});
semilogy(setting.SNR_dB,sqrt(meanDTFT.  MAP_MSE_Signal),spec.DTFT_MAP{:});
semilogy(setting.SNR_dB,sqrt(meanMUSIC.trueK_MSE_Signal),spec.MUSIC_trueK{:});
semilogy(setting.SNR_dB,sqrt(meanMUSIC.  AIC_MSE_Signal),spec.MUSIC_AIC{:});
semilogy(setting.SNR_dB,sqrt(meanMUSIC.  MAP_MSE_Signal),spec.MUSIC_MAP{:});
%-----------------------------------------
xlabel('SNR (dB)');
ylabel('RMSE of ML estimate of amplitudes');

set(gca, 'YScale', 'log')

legend
movegui(h.ML,'northwest')
%//////////////////////////////////////////////////////////////////////////

h.MAP = figure;
hold on;
%-----------------------------------------
semilogy(setting.SNR_dB,sqrt(meanDTFT.trueK_MSE_SignalTau),spec.DTFT_trueK{:});
semilogy(setting.SNR_dB,sqrt(meanDTFT.  AIC_MSE_SignalTau),spec.DTFT_AIC{:});
semilogy(setting.SNR_dB,sqrt(meanDTFT.  MAP_MSE_SignalTau),spec.DTFT_MAP{:});
semilogy(setting.SNR_dB,sqrt(meanMUSIC.trueK_MSE_SignalTau),spec.MUSIC_trueK{:});
semilogy(setting.SNR_dB,sqrt(meanMUSIC.  AIC_MSE_SignalTau),spec.MUSIC_AIC{:});
semilogy(setting.SNR_dB,sqrt(meanMUSIC.  MAP_MSE_SignalTau),spec.MUSIC_MAP{:});
%-----------------------------------------
xlabel('SNR (dB)');
ylabel('RMSE of MAP estimate of amplitudes');

set(gca, 'YScale', 'log')

legend
movegui(h.MAP,'southwest')

%//////////////////////////////////////////////////////////////////////////
