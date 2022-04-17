% ECG_ABP_PPG_RRFinder - takes an input ECG signal, PPG signal and ABP signal, calculates wavelet sequence,
% fuses wavelets, detects peaks then returns RR interval sequence
% Author - Shane Malone, Student Number 17308336
% Credit - Arlene John for algorithm and basic code structure

function [RRIntervalSet, peakLocations, final_sig_tab, signal2plot, ecgsqi_orig, ppgsqi_orig, abpsqi_orig, ppg_wav_win, ecg_wav_win, abp_wav_win, ppgfin_sig_tab, ecgfin_sig_tab, abpfin_sig_tab]= ECG_ABP_PPG_RRFinder(ecgInput,ppgInput,abpInput,Fs,ppgDelay, abpDelay)
signal(:,1)=double(ecgInput(round(2*Fs):end-round(2*Fs)+1));
signal(:,1)=signal(:,1)-mean(signal(:,1));

signal(:,2)=double(ppgInput(round(2*Fs)-ppgDelay:end-round(2*Fs)-ppgDelay+1));

signal(:,3)=double(abpInput(round(2*Fs)-abpDelay:end-round(2*Fs)-abpDelay+1));

% load template curves
load('PPGData.mat')
load('ECGData.mat');
load('ABP.mat')

% ECG SQI
ecgsignal=signal(:,1)';

ecgCurve=double(ECG.signal);
curvefreq=ECG.Fs;
ecgsqi=csqi_calc(ecgsignal,Fs,ecgCurve,curvefreq); 

% PPG SQI
ppgsignal=signal(:,2);
ppgCurve=double(PPG.signal);
curvefreq=PPG.Fs;
ppgsqi=pcsqi_calc(ppgsignal,Fs,ppgCurve,curvefreq);
ppgsignal=ppgsignal';

% ABP SQI
abpsignal=signal(:,3);
abpCurve=double(ABP.signal);
curvefreq=ABP.Fs;
abpsqi = Acsqi_calc(abpsignal,Fs,abpCurve,curvefreq);
abpsignal=abpsignal';

% Replace NaN with 0
if sum(isnan(ecgsqi)) > 0
    ecgsqi = zeros(length(ecgsqi), 1);
end
if sum(isnan(ppgsqi)) > 0
    ppgsqi = zeros(length(ppgsqi), 1);
end
if sum(isnan(abpsqi)) > 0
    abpsqi = zeros(length(abpsqi), 1);
end

ppgsqi_orig=ppgsqi;
ecgsqi_orig=ecgsqi;
abpsqi_orig=abpsqi;

% Normalizing the SQIs to 1
sum_sqi=ecgsqi+ppgsqi+abpsqi;
ecgsqi(sum_sqi==0)=0.34;
ppgsqi(sum_sqi==0)=0.33;
abpsqi(sum_sqi==0)=0.33;
sum_sqi=ecgsqi+ppgsqi+abpsqi;

ecgsqi=ecgsqi./sum_sqi;
ppgsqi=ppgsqi./sum_sqi;
abpsqi=abpsqi./sum_sqi;

ecgsqi=ecgsqi';
ppgsqi=ppgsqi';
abpsqi=abpsqi';


window_size=round(3*Fs)+1; % choosing window size for processing
k=1;

ann=[0 0]; % pre allocating annotation variable for speed

%  preparing the windowed wavelet coefficients
ppg_wav_win = zeros(1, 1);
ecg_wav_win = zeros(1, 1);
abp_wav_win = zeros(1, 1);

final_sig_tab=zeros(length(ecgsignal),1);
ppgfin_sig_tab=zeros(length(ecgsignal),1);
ecgfin_sig_tab=zeros(length(ecgsignal),1);
abpfin_sig_tab=zeros(length(ecgsignal),1);


while k<=length(ecgsignal)
    ecgwindow=ecgsignal(k:min(k+window_size-1,length(ecgsignal)));
    len=length(ecgwindow);
    
    % If window size lower than selected window size
    if length(ecgwindow)<window_size/2
        k=k-(window_size-length(ecgwindow));
        ecgwindow=ecgsignal(k:min(k+window_size-1,length(ecgsignal)));
        len=length(ecgwindow);
    end
    
    % Select ecgwindow, ppgwindow, abpwindow, sqi windows and append 0s in the beginning and prepend 0s at the end
    ecgwindow=[zeros(1,Fs) ecgwindow zeros(1,Fs)];
    ppgwindow=[zeros(1,Fs) ppgsignal(k:min(k+window_size-1,length(ppgsignal))) zeros(1,Fs)];
    abpwindow=[zeros(1,Fs) abpsignal(k:min(k+window_size-1,length(abpsignal))) zeros(1,Fs)];
    ecgsqiwindow=[zeros(1,Fs) ecgsqi(k:min(k+window_size-1,length(ecgsignal))) zeros(1,Fs)];
    ppgsqiwindow=[zeros(1,Fs) ppgsqi(k:min(k+window_size-1,length(ppgsignal))) zeros(1,Fs)];
    abpsqiwindow=[zeros(1,Fs) abpsqi(k:min(k+window_size-1,length(abpsignal))) zeros(1,Fs)];
   
    ppg_sig=[0];
    abp_sig=[0];
    ecg_sig=[0];
         
    if mod(length(ecgwindow),8)~=0 % adjusting length for wavelet transform
        ecgwindow=[zeros(1,floor((8-mod(length(ecgwindow),8))/2)) ecgwindow zeros(1,ceil((8-mod(length(ecgwindow),8))/2))];
        ecgsqiwindow=[zeros(1,floor((8-mod(length(ecgsqiwindow),8))/2)) ecgsqiwindow zeros(1,ceil((8-mod(length(ecgsqiwindow),8))/2))];
    end
    if mod(length(ppgwindow),8)~=0
        ppgwindow=[zeros(1,floor((8-mod(length(ppgwindow),8))/2)) ppgwindow zeros(1,ceil((8-mod(length(ppgwindow),8))/2))];
        ppgsqiwindow=[zeros(1,floor((8-mod(length(ppgsqiwindow),8))/2)) ppgsqiwindow zeros(1,ceil((8-mod(length(ppgsqiwindow),8))/2))];
    end
    if mod(length(abpwindow),8)~=0
        abpwindow=[zeros(1,floor((8-mod(length(abpwindow),8))/2)) abpwindow zeros(1,ceil((8-mod(length(abpwindow),8))/2))];
        abpsqiwindow=[zeros(1,floor((8-mod(length(abpsqiwindow),8))/2)) abpsqiwindow zeros(1,ceil((8-mod(length(abpsqiwindow),8))/2))];
    end

    % ECG Wavelet transform
    [~,swd] = swt(ecgwindow,3,'bior6.8');
    ecg_sig=swd(2,:)+swd(3,:);
    ecg_sig=ecg_sig(round((length(ecg_sig)-len)/2+1):round((length(ecg_sig)+len)/2)); % adjust length to window length
    ecgsqiwindow=ecgsqiwindow(round((length(ecgsqiwindow)-len)/2+1):round((length(ecgsqiwindow)+len)/2));

    % PPG Wavelet transform
    [~,swd] = swt(ppgwindow,3,'sym8');
    ppg_sig=swd(2,:)+swd(3,:);
    if sum(isnan(ppg_sig)) > 0
        ppg_sig = zeros(1, length(ppg_sig));
    end
    ppg_sig=ppg_sig(round((length(ppg_sig)-len)/2+1):round((length(ppg_sig)+len)/2));
    ppgsqiwindow=ppgsqiwindow(round((length(ppgsqiwindow)-len)/2+1):round((length(ppgsqiwindow)+len)/2));

    % ABP Wavelet transform
    [~,swd] = swt(abpwindow,3,'bior6.8');
    abp_sig=swd(2,:)+swd(3,:);
    if sum(isnan(abp_sig)) > 0
        abp_sig = zeros(1, length(abp_sig));
    end
    abp_sig=abp_sig(round((length(abp_sig)-len)/2+1):round((length(abp_sig)+len)/2));
    abpsqiwindow=abpsqiwindow(round((length(abpsqiwindow)-len)/2+1):round((length(abpsqiwindow)+len)/2));
    
    % ADDITIVE FUSION
    abpScale = 1/20;
    ppgfin_sig = ppg_sig.*ppgsqiwindow;
    ecgfin_sig = ecg_sig.*ecgsqiwindow;
    abpfin_sig = abpScale*abp_sig.*abpsqiwindow;

    ppg_wav_win=[ppg_wav_win ppg_sig]; % variables to track wavelet sequence
    abp_wav_win=[abpScale*abp_wav_win abp_sig];
    ecg_wav_win=[ecg_wav_win ecg_sig];
    
    final_sig = ppgfin_sig+ecgfin_sig+abpfin_sig;

    ppgfin_sig_tab(k:k+length(final_sig)-1)=ppgfin_sig;
    abpfin_sig_tab(k:k+length(final_sig)-1)=abpfin_sig;
    ecgfin_sig_tab(k:k+length(final_sig)-1)=ecgfin_sig;
    final_sig_tab(k:k+length(final_sig)-1)=final_sig; % Tab of signal used for peak detection
       
% peak detection
    final_sig(1:10) = 0; % This is to ignore the peak in wavelet sequence at start of signal
    [pksmax,locsmax] = findpeaks(final_sig,'MINPEAKDISTANCE',round(0.3*Fs)); % find peaks every 0.3 s
    Vmax=0.5*mean(pksmax); % set a threshold based on the peaks
    maximatokeep=pksmax>=Vmax; % only peaks above the threshold
    locsmaxima=locsmax([maximatokeep]); % locations of selected peaks
    maxmin_locs=sort(locsmaxima); % sort the locations
    maxmin_val=final_sig(maxmin_locs); % find the peaks at the locations 
    l=1;
    clear hbloc
    for i=1:length(maxmin_val)
        if(maxmin_val(i)>=Vmax)
            hbloc(l)=k+maxmin_locs(i)-1; % correct the peak annotation value within this window to that of the signal as a whole
            l=l+1;
        end
    end
    if l>1
        ann=[ann hbloc];
    end
    
    k=k+window_size;
end

% extract the detected peaks
ann=ann(3:end);
ann = unique(ann);
peakLocations=[0]; % pre-allocate secondary annotation locations

for k=1:length(ann)-2
    if peakLocations(end)>=ann(k)% in case the position at the last position of sec locs is greater than that of ann
        continue
    end
    rrinterval=(ann(k+1)-ann(k))/Fs; % calculate rrinterval
    t=0;
    if rrinterval > 1.2 % if rrinterval is greter than 1.2 s, means missed beat
        signalwindow=final_sig_tab(ann(k)+round(rrinterval*Fs/3):ann(k+1)-round(rrinterval*Fs/3));
        if sum(signalwindow)~=0
            [~,idx_val]=max(signalwindow);
            hbloc=ann(k)+idx_val+round(rrinterval*Fs/3);
            peakLocations=[peakLocations ann(k) hbloc];
            t=1;
        end  
    end
    if rrinterval<0.3 % if rrinterval is less than 0.4s,means false beat
        if rrinterval<0.15 % if two peaks within 0.15 s, take average
            hbloc=ann(k+1);
            peakLocations=[peakLocations hbloc];
            t=1;
        else
            rrinterval2=(ann(k+2)-ann(k))/Fs; % check interval between the next to next beat
            if rrinterval2 < 1.2 % if less than 1.2, omit k+1 location
                peakLocations=[peakLocations ann(k) ann(k+2)];  
                t=1;
            else
                peakLocations=[peakLocations ann(k)];% else keep it for next round
                t=1;
            end
        end
    end
    if t==0 % if no corrections made, just add at the end
       peakLocations=[peakLocations ann(k)]; 
    end
end

peakLocations=[peakLocations ann(k+1)];
peakLocations=peakLocations(2:end);
peakLocations=unique(peakLocations); % final annotations
ppg_wav_win=ppg_wav_win(2:end);
abp_wav_win=abp_wav_win(2:end);
ecg_wav_win=ecg_wav_win(2:end);
signal2plot=[ecgsignal;ppgsignal;abpsignal]'; % for plotting

% Finally get rr interval
RRIntervalSet = zeros(1, length(peakLocations)-1);
for k = 1:length(peakLocations)-1
    rrinterval = (peakLocations(k+1) - peakLocations(k))/Fs; % calculate rrinterval
    RRIntervalSet(k) = rrinterval;
end