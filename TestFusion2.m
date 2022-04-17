SignalLocs = readtable('ECG_PPG_SignalLocations.csv');
row = 42;
Fs = 125;

RecordName = SignalLocs{row, 'Record'};
ECGLoc = SignalLocs{row, 'ECG'};
PPGLoc = SignalLocs{row, 'PPG'};
ppgDelay = SignalLocs{row, 'Delay'};
%ppgDelay = 50;

startTime = 60000;
[signal,~,~]=rdsamp(RecordName{1},[],(startTime+10)*Fs,startTime*Fs,1); %%%%% read digital signal values

ecgInput = signal(:,ECGLoc);
ppgInput = signal(:,PPGLoc);
%ecgInput = awgn(ecgInput, -10, 'measured');
%ppgInput = awgn(ppgInput, -10, 'measured');


[RRIntervalSet, sec_locs, final_sig_tab, signal2plot, ecgsqi_orig, ppgsqi_orig, ppg_wav_win, ecg_wav_win, ppgfin_sig_tab, ecgfin_sig_tab] = ECG_PPG_RRFinder(ecgInput,ppgInput,Fs,ppgDelay);

% plotting
time_base=1:length(signal2plot);
time_base=time_base/Fs;
%subplot = @(m,n,p) subtightplot(m,n,p,[0.05,0.01],[0.05,0.05],[0.1,0.03]);
figure
ax1=subplot(9,1,1);
plot(time_base,signal2plot(:,1),'r','Linewidth',1.2)
title('ECG signal')

ax2=subplot(9,1,2);
plot(time_base,signal2plot(:,2),'b','Linewidth',1.2)
title('PPG signal')

ax3=subplot(9,1,3);
plot(time_base,ecgsqi_orig,'r','Linewidth',1.2)
title('ECG lincSQI')

ax4=subplot(9,1,4);
plot(time_base,ppgsqi_orig,'b','Linewidth',1.2)
title('PPG lincSQI')

wav_time_base = (1:length(ecg_wav_win))/Fs;
ax5=subplot(9,1,5);
plot(wav_time_base,ecg_wav_win,'r','Linewidth',1.2)
title('ECG d2+d3 feature sequence')

ax6=subplot(9,1,6);
plot(wav_time_base,ppg_wav_win,'b','Linewidth',1.2)
title('PPG d2+d3 feature sequence')

ax7=subplot(9,1,7);
plot(time_base,ecgfin_sig_tab,'r','Linewidth',1.2)
title('ECG lincSQI * ECG d2+d3 feature sequence')

ax8=subplot(9,1,8);
plot(time_base,ppgfin_sig_tab,'b','Linewidth',1.2)
title('PPG lincSQI* PPG d2+d3 feature sequence')

ax9=subplot(9,1,9);
plot(time_base,final_sig_tab,'k','Linewidth',1.2)
hold on
plot(sec_locs./125,zeros(length(sec_locs)),'mo','Linewidth',1.2)

title('Fused sequence with location from fusion algorithm')

p1=get(ax4,'position');
p2=get(ax5,'position');
height=p1(2)+p1(4)-p2(2);
h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
h_label=ylabel('Amplitude (a.u.)','visible','on');