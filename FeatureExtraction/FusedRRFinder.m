function [RRIntervalSet, sec_locs, final_sig_tab, signal2plot, ecgsqi_orig, ppgsqi_orig, ppg_wav_win, ecg_wav_win, ppgfin_sig_tab, ecgfin_sig_tab]= FusedRRFinder(ecgInput,ppgInput,Fs,delay)

%%%%signalorig -signal
%%%%Fs- sampling frequnecy of signal
signal(:,1)=double(ecgInput(round(2*Fs):end-round(2*Fs)+1)); %%%% 1st column - ECG signal ; start after 2 seconds
signal(:,1)=signal(:,1)-mean(signal(:,1));

signal(:,2)=double(ppgInput(round(2*Fs)-delay:end-round(2*Fs)-delay+1));

ecgsignal=signal(:,1)';

%%%%%load template curves
load('ppgdata1.mat');
load('ecgdata1.mat');

%%%%extracting the curve template
ecgcurve=double(ECG.signal);
curvefreq=ECG.Fs;
%%%%calculating ECG SQI
ecgsqi=csqi_calc(ecgsignal,Fs,ecgcurve,curvefreq); 

    

    
%%%%PPGSIGNAL Channel 2

ppgsignal=signal(:,2)-mean(signal(:,2)); % Takes away mean to have centred on zero?


%%%%%extracting curve template

curve=double(PPG.signal);
curvefreq=PPG.Fs;

%%%%calculating PPG SQI
ppgsqi=pcsqi_calc(ppgsignal,Fs,curve,curvefreq);
ppgsqi_orig=ppgsqi;
ecgsqi_orig=ecgsqi;

    
%%%%%Normalizing the SQIs to 1
sum_sqi=ecgsqi+ppgsqi;
ecgsqi(sum_sqi==0)=0.5;
ppgsqi(sum_sqi==0)=0.5;
sum_sqi=ecgsqi+ppgsqi;
ecgsqi=ecgsqi./sum_sqi;
ppgsqi=ppgsqi./sum_sqi;
ecgsqi=ecgsqi';
ppgsqi=ppgsqi';
ppgsignal=ppgsignal';



window_size=round(3*Fs)+1; %%%%choosing window size for processing
k=1;

ann=[0 0]; %%%pre allocating annotation variable fo speed

%%%%% preparing the windowed wavelet coefficients (this portion of the code will be different for
%%%%%%%%%%%respiration rate)
ppg_wav_win=[0];
ecg_wav_win=[0];
final_sig_tab=zeros(length(ecgsignal),1); %%%%tab of the final processed signal
ppgfin_sig_tab=zeros(length(ecgsignal),1);
ecgfin_sig_tab=zeros(length(ecgsignal),1);
while k<=length(ecgsignal)
    ecgwindow=ecgsignal(k:min(k+window_size-1,length(ecgsignal))); %%%%ecg window
    len=length(ecgwindow);
    
    %%%%if window size lower than selected window size
    if length(ecgwindow)<window_size/2
        k=k-(window_size-length(ecgwindow));
        ecgwindow=ecgsignal(k:min(k+window_size-1,length(ecgsignal)));
        len=length(ecgwindow);
    end
    
    %%%%%select ecgwindow, ppgwindow, sqi windows and append 0s in the
    %%%%%beginning and prepend 0s at the end
    ecgwindow=[zeros(1,Fs) ecgwindow zeros(1,Fs)];
    ppgwindow=[zeros(1,Fs) ppgsignal(k:min(k+window_size-1,length(ppgsignal))) zeros(1,Fs)];
    ecgsqiwindow=[zeros(1,Fs) ecgsqi(k:min(k+window_size-1,length(ecgsignal))) zeros(1,Fs)];
    ppgsqiwindow=[zeros(1,Fs) ppgsqi(k:min(k+window_size-1,length(ppgsignal))) zeros(1,Fs)];
   
    ppg_sig=[0];
    ecg_sig=[0];
         
    if mod(length(ecgwindow),8)~=0 %%%%adjusting length for wavelet transform
        ecgwindow=[zeros(1,floor((8-mod(length(ecgwindow),8))/2)) ecgwindow zeros(1,ceil((8-mod(length(ecgwindow),8))/2))];
        ecgsqiwindow=[zeros(1,floor((8-mod(length(ecgsqiwindow),8))/2)) ecgsqiwindow zeros(1,ceil((8-mod(length(ecgsqiwindow),8))/2))];
    end
    if mod(length(ppgwindow),8)~=0
        ppgwindow=[zeros(1,floor((8-mod(length(ppgwindow),8))/2)) ppgwindow zeros(1,ceil((8-mod(length(ppgwindow),8))/2))];
        ppgsqiwindow=[zeros(1,floor((8-mod(length(ppgsqiwindow),8))/2)) ppgsqiwindow zeros(1,ceil((8-mod(length(ppgsqiwindow),8))/2))];
    end
    %%%%%%%%%%%%%%ECG Wavelet transform
    [~,swd] = swt(ecgwindow,3,'bior6.8');  %%%%stationary wavelet tranform
%             ecg_sig=abs(swd(2,:)).*abs(swd(3,:));
%     ecg_sig=abs(swd(2,:))+abs(swd(3,:));
    ecg_sig=swd(2,:)+swd(3,:);
    ecg_sig=ecg_sig(round((length(ecg_sig)-len)/2+1):round((length(ecg_sig)+len)/2)); %%%%adjust length to window length
    ecgsqiwindow=ecgsqiwindow(round((length(ecgsqiwindow)-len)/2+1):round((length(ecgsqiwindow)+len)/2));

    %%%%%%%%%%%%%%PPG Wavelet transform
    [~,swd] = swt(ppgwindow,3,'sym8');
%   ppg_sig=12.4545.*abs(swd(2,:)).*abs(swd(3,:));
%     ppg_sig=(abs(swd(2,:))+abs(swd(3,:)));
    ppg_sig=swd(2,:)+swd(3,:);
    ppg_sig=ppg_sig(round((length(ppg_sig)-len)/2+1):round((length(ppg_sig)+len)/2));ppgsqiwindow=ppgsqiwindow(round((length(ppgsqiwindow)-len)/2+1):round((length(ppgsqiwindow)+len)/2));
    
%%%%%%%ADDITIVE FUSION
  
   ppgfin_sig=ppg_sig.*4.*ppgsqiwindow; %%ppg wavelet sequence amplitude is almost 4 times smaller than ecg
   ecgfin_sig=ecg_sig.*ecgsqiwindow;
   ppg_wav_win=[ppg_wav_win ppg_sig.*4]; %%variables to track wavelet sequence
   ecg_wav_win=[ecg_wav_win ecg_sig];
    
    final_sig=ppgfin_sig+ecgfin_sig;

    ppgfin_sig_tab(k:k+length(final_sig)-1)=ppgfin_sig;
    ecgfin_sig_tab(k:k+length(final_sig)-1)=ecgfin_sig;
    final_sig_tab(k:k+length(final_sig)-1)=final_sig; %%%%Tab of signal used for peak detection
       
%%%%%%%%%%%peak detection (this portion of the code will be different for
%%%%%%%%%%%respiration rate)
    [pksmax,locsmax] = findpeaks(final_sig,'MINPEAKDISTANCE',round(0.3*Fs)); %%%%find peaks every 0.3 s
    Vmax=0.5*mean(pksmax); %%%%set a threshold based on the peaks
    maximatokeep=pksmax>=Vmax; %%%%%only peaks above the threshold
    locsmaxima=locsmax([maximatokeep]); %%%%%locations of selected peaks
    maxmin_locs=sort(locsmaxima); %%%%sort the locations
    maxmin_val=final_sig(maxmin_locs); %%%%find the peaks at the locations 
    l=1;
    clear hbloc
    for i=1:length(maxmin_val)
        if(maxmin_val(i)>=Vmax)
            hbloc(l)=k+maxmin_locs(i)-1; %%%%correct the peak annotation value within this window to that of the signal as a whole
            l=l+1;
        end
    end
    if l>1
        ann=[ann hbloc];
    end
    
    k=k+window_size;
end

%%%%extract the detected peaks
ann=ann(3:end);
%%%%check for false negatives and positives and correct%%%%%%

sec_locs=[0]; %%%% pre-allocate secondary annotation locations

for k=1:length(ann)-2
    if sec_locs(end)>=ann(k)%%%%in case the position at the last position of sec locs is greater than that of ann
        continue
    end
    rrinterval=(ann(k+1)-ann(k))/Fs; %%%%calculate rrinterval
    t=0;
    if rrinterval > 1.2 %%%if rrinterval is greter than 1.2 s, means missed beat
        signalwindow=final_sig_tab(ann(k)+round(rrinterval*Fs/3):ann(k+1)-round(rrinterval*Fs/3));
        if sum(signalwindow)~=0
            [~,idx_val]=max(signalwindow);
            hbloc=ann(k)+idx_val+round(rrinterval*Fs/3);
            sec_locs=[sec_locs ann(k) hbloc];
            t=1;
        end  
    end
    if rrinterval<0.3 %%%%if rrinterval is less than 0.4s,means false beat
        if rrinterval<0.1 %%%%if two peaks within 0.1 s, take average
            hbloc=ann(k+1);
            sec_locs=[sec_locs hbloc];
            t=1;
        else
            rrinterval2=(ann(k+2)-ann(k))/Fs; %%%%check interval between the next to next beat
            if rrinterval2 < 1.2 %%%%if less than 1.2, omit k+1 location
                sec_locs=[sec_locs ann(k) ann(k+2)];  
                t=1;
            else
                sec_locs=[sec_locs ann(k)];%%%%else keep it for next round
                t=1;
            end
        end
    end
    if t==0 %%%if no corrections made, just add at the end
       sec_locs=[sec_locs ann(k)]; 
    end
end

sec_locs=[sec_locs ann(k+1)];
sec_locs=sec_locs(2:end);
sec_locs=unique(sec_locs); %%%%final annotations
ppg_wav_win=ppg_wav_win(2:end);
ecg_wav_win=ecg_wav_win(2:end);
signal2plot=[ecgsignal;ppgsignal]'; %%%%for plotting purpo

% Finally get rr interval
RRIntervalSet = zeros(1, length(sec_locs)-1);
for k = 1:length(sec_locs)-1
    rrinterval = (sec_locs(k+1) - sec_locs(k))/Fs; %%%%calculate rrinterval
    RRIntervalSet(k) = rrinterval;
end