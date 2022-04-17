SignalLocs = readtable('ECG_PPG_ABP_SignalLocationsUpdated.csv');
row = 1;
Fs = 125;

RecordName = SignalLocs{row, 'Record'};
ECGLoc = SignalLocs{row, 'ECG'};
PPGLoc = SignalLocs{row, 'PPG'};
ABPLoc = SignalLocs{row, 'ABP'};
ppgDelay = SignalLocs{row, 'Delay'};
abpDelay = SignalLocs{row, 'ABPDelay'};
%abpDelay = -30;

startTime = 60000;
[signal,~,~]=rdsamp(RecordName{1},[],(startTime+15)*Fs,startTime*Fs,1); % read digital signal values

Length = 10*Fs;
SNR = -10;
ecgStart = floor(2/8*Length);
ecgEnd = floor(3/8*Length);
ppgStart = floor(3/8*Length);
ppgEnd = floor(4/8*Length);
abpStart = floor(4/8*Length);
abpEnd = floor(5/8*Length);

ecgInput = signal(:,ECGLoc);
% ecgInput = ecgInput - min(ecgInput(:)); % Normalise ppg
% ecgInput = ecgInput ./ max(ecgInput(:));

ppgInput = signal(:,PPGLoc);
ppgInput = ppgInput - min(ppgInput(:)); % Normalise ppg
ppgInput = ppgInput ./ max(ppgInput(:));

abpInput = signal(:,ABPLoc);


ecgInput(ecgStart:ecgEnd) = awgn(ecgInput(ecgStart:ecgEnd), SNR, 'measured');
ppgInput(ppgStart:ppgEnd) = awgn(ppgInput(ppgStart:ppgEnd), SNR, 'measured');
abpInput(abpStart:abpEnd) = awgn(abpInput(abpStart:abpEnd), SNR, 'measured');



[RRIntervalSet, sec_locs, final_sig_tab, signal2plot, ecgsqi_orig, ppgsqi_orig, abpsqi_orig, ppg_wav_win, ecg_wav_win, abp_wav_win, ppgfin_sig_tab, ecgfin_sig_tab, abpfin_sig_tab] = ECG_ABP_PPG_RRFinder(ecgInput,ppgInput,abpInput,Fs,ppgDelay,abpDelay);

% plotting
time_base=1:length(signal2plot);
time_base=time_base/Fs;
%subplot = @(m,n,p) subtightplot(m,n,p,[0.05,0.01],[0.05,0.05],[0.1,0.03]);
figure
ax1=subplot(10,1,1);
plot(time_base,signal2plot(:,1),'r','Linewidth',1.2)
title('ECG signal')

ax2=subplot(10,1,2);
plot(time_base,signal2plot(:,2),'b','Linewidth',1.2)
title('PPG signal')

ax3=subplot(10,1,3);
plot(time_base,signal2plot(:,3),'b','Linewidth',1.2)
title('ABP signal')

ax4=subplot(10,1,4);
plot(time_base,ecgsqi_orig,'r','Linewidth',1.2)
title('ECG lincSQI')
ylim([0, 20])

ax5=subplot(10,1,5);
plot(time_base,ppgsqi_orig,'b','Linewidth',1.2)
title('PPG lincSQI')
ylim([0, 20])

ax6=subplot(10,1,6);
plot(time_base,abpsqi_orig,'b','Linewidth',1.2)
title('ABP lincSQI')
ylim([0, 20])

ax7=subplot(10,1,7);
plot(time_base,ecgfin_sig_tab,'r','Linewidth',1.2)
title('ECG lincSQI * ECG d2+d3 feature sequence')
ylim([-0.5, 0.5])

ax8=subplot(10,1,8);
plot(time_base,ppgfin_sig_tab,'b','Linewidth',1.2)
title('PPG lincSQI* PPG d2+d3 feature sequence')
ylim([-0.5, 0.5])

ax9=subplot(10,1,9);
plot(time_base,abpfin_sig_tab,'b','Linewidth',1.2)
title('ABP lincSQI* ABP d2+d3 feature sequence')
ylim([-0.5, 0.5])

ax10=subplot(10,1,10);
plot(time_base,final_sig_tab,'k','Linewidth',1.2)
hold on
plot(sec_locs./125,zeros(length(sec_locs)),'mo','Linewidth',1.2)

title('Fused sequence with location from fusion algorithm')

p1=get(ax4,'position');
p2=get(ax5,'position');
height=p1(2)+p1(4)-p2(2);
h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
h_label=ylabel('Amplitude (a.u.)','visible','on');


% % plotting
% time_base=1:length(signal2plot);
% time_base=time_base/Fs;
% %subplot = @(m,n,p) subtightplot(m,n,p,[0.05,0.01],[0.05,0.05],[0.1,0.03]);
% figure
% ax1=subplot(13,1,1);
% plot(time_base,signal2plot(:,1),'r','Linewidth',1.2)
% title('ECG signal')
% 
% ax2=subplot(13,1,2);
% plot(time_base,signal2plot(:,2),'b','Linewidth',1.2)
% title('PPG signal')
% 
% ax3=subplot(13,1,3);
% plot(time_base,signal2plot(:,3),'b','Linewidth',1.2)
% title('ABP signal')
% 
% ax4=subplot(13,1,4);
% plot(time_base,ecgsqi_orig,'r','Linewidth',1.2)
% title('ECG lincSQI')
% 
% ax5=subplot(13,1,5);
% plot(time_base,ppgsqi_orig,'b','Linewidth',1.2)
% title('PPG lincSQI')
% 
% ax6=subplot(13,1,6);
% plot(time_base,abpsqi_orig,'b','Linewidth',1.2)
% title('ABP lincSQI')
% 
% wav_time_base = (1:length(ecg_wav_win))/Fs;
% ax7=subplot(13,1,7);
% plot(wav_time_base,ecg_wav_win,'r','Linewidth',1.2)
% title('ECG d2+d3 feature sequence')
% 
% ax8=subplot(13,1,8);
% plot(wav_time_base,ppg_wav_win,'b','Linewidth',1.2)
% title('PPG d2+d3 feature sequence')
% 
% ax9=subplot(13,1,9);
% plot(wav_time_base,abp_wav_win,'b','Linewidth',1.2)
% title('ABP d2+d3 feature sequence')
% 
% ax10=subplot(13,1,10);
% plot(time_base,ecgfin_sig_tab,'r','Linewidth',1.2)
% title('ECG lincSQI * ECG d2+d3 feature sequence')
% 
% ax11=subplot(13,1,11);
% plot(time_base,ppgfin_sig_tab,'b','Linewidth',1.2)
% title('PPG lincSQI* PPG d2+d3 feature sequence')
% 
% ax12=subplot(13,1,12);
% plot(time_base,abpfin_sig_tab,'b','Linewidth',1.2)
% title('ABP lincSQI* ABP d2+d3 feature sequence')
% 
% ax13=subplot(13,1,13);
% plot(time_base,final_sig_tab,'k','Linewidth',1.2)
% hold on
% plot(sec_locs./125,zeros(length(sec_locs)),'mo','Linewidth',1.2)
% 
% title('Fused sequence with location from fusion algorithm')
% 
% p1=get(ax4,'position');
% p2=get(ax5,'position');
% height=p1(2)+p1(4)-p2(2);
% h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
% h_label=ylabel('Amplitude (a.u.)','visible','on');