% PPGRRFinder - takes an input PPG signal, calculates wavelet sequence,
% detects peaks then returns RR interval sequence
% Author - Shane Malone, Student Number 17308336
% Credit - Arlene John for algorithm and basic code structure

function [RRIntervalSet, peakLocations, final_sig_tab, signal2plot, ppg_wav_win, ppgfin_sig_tab]= PPGRRFinder(ppgInput,Fs)

signal(:,1)=double(ppgInput(round(2*Fs):end-round(2*Fs)+1));

ppgsignal=signal(:,1)-mean(signal(:,1)); % Takes away mean to have centred on zero

ppgsignal=ppgsignal';



window_size=round(3*Fs)+1; % choosing window size for processing
k=1;

ann=[0 0]; % pre allocating annotation variable fo speed


ppg_wav_win=[0];

final_sig_tab=zeros(length(ppgsignal),1); % tab of the final processed signal
ppgfin_sig_tab=zeros(length(ppgsignal),1);

while k<=length(ppgsignal)
    ecgwindow=ppgsignal(k:min(k+window_size-1,length(ppgsignal))); % ecg window
    len=length(ecgwindow);
    
    % if window size lower than selected window size
    if length(ecgwindow)<window_size/2
        k=k-(window_size-length(ecgwindow));
        ecgwindow=ppgsignal(k:min(k+window_size-1,length(ppgsignal)));
        len=length(ecgwindow);
    end
    
    % select ecgwindow, ppgwindow, sqi windows and append 0s in the
    % beginning and prepend 0s at the end
    ecgwindow=[zeros(1,Fs) ecgwindow zeros(1,Fs)];
    ppgwindow=[zeros(1,Fs) ppgsignal(k:min(k+window_size-1,length(ppgsignal))) zeros(1,Fs)];
   
    ppg_sig=[0];
         
    if mod(length(ecgwindow),8)~=0 % adjusting length for wavelet transform
        ecgwindow=[zeros(1,floor((8-mod(length(ecgwindow),8))/2)) ecgwindow zeros(1,ceil((8-mod(length(ecgwindow),8))/2))];
    end
    if mod(length(ppgwindow),8)~=0
        ppgwindow=[zeros(1,floor((8-mod(length(ppgwindow),8))/2)) ppgwindow zeros(1,ceil((8-mod(length(ppgwindow),8))/2))];
    end

    % PPG Wavelet transform
    [~,swd] = swt(ppgwindow,3,'sym8');
    ppg_sig=swd(2,:)+swd(3,:);
    ppg_sig=ppg_sig(round((length(ppg_sig)-len)/2+1):round((length(ppg_sig)+len)/2));
    
    % ADDITIVE FUSION
  
    ppgfin_sig=ppg_sig.*4; % scale ppg to match ECG amplitude
    ppg_wav_win=[ppg_wav_win ppg_sig.*4]; % variables to track wavelet sequence
    
    final_sig=ppgfin_sig;

    ppgfin_sig_tab(k:k+length(final_sig)-1)=ppgfin_sig;
    final_sig_tab(k:k+length(final_sig)-1)=final_sig; % Tab of signal used for peak detection
       
    % Peak Detection
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
            hbloc(l)=k+maxmin_locs(i)-1;
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
    if rrinterval > 1.2 % if rrinterval is greater than 1.2 s, means missed beat
        signalwindow=final_sig_tab(ann(k)+round(rrinterval*Fs/3):ann(k+1)-round(rrinterval*Fs/3));
        if sum(signalwindow)~=0
            [~,idx_val]=max(signalwindow);
            hbloc=ann(k)+idx_val+round(rrinterval*Fs/3);
            peakLocations=[peakLocations ann(k) hbloc];
            t=1;
        end  
    end
    if rrinterval<0.3 % if rrinterval is less than 0.3s,means false beat
        if rrinterval<0.2 % if two peaks within 0.1 s, take average
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
ppg_wav_win=ppg_wav_win(2:end);
signal2plot=[ppgsignal]'; % for plotting purposes

% Finally get rr interval
RRIntervalSet = zeros(1, length(peakLocations)-1);
for k = 1:length(peakLocations)-1
    rrinterval = (peakLocations(k+1) - peakLocations(k))/Fs; % calculate rrinterval
    RRIntervalSet(k) = rrinterval;
end