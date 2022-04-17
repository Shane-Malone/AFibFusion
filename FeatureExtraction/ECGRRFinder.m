% ECGRRFinder - takes an input ECG signal, calculates wavelet sequence,
% detects peaks then returns RR interval sequence
% Author - Shane Malone, Student Number 17308336
% Credit - Arlene John for algorithm and basic code structure

function [RRIntervalSet, peakLocations, finalSignalTab, signal2plot, ecgWavWindow, ecgSignalTab]= ECGRRFinder(signalOrig, Fs)


signal(:,1) = double(signalOrig(round(2*Fs):end-round(2*Fs)+1));
signal(:,1) = signal(:,1)-mean(signal(:,1));

ecgSignal = signal(:,1)';
lengthInput = length(ecgSignal);

% Load template curves
load('ecgdata1.mat');

window_size = round(3*Fs)+1; % Choosing window size for processing
k = 1;

annotation = [0 0]; % pre allocating annotation variable for speed

% Preparing the windowed wavelet coefficients
ecgWavWindow = [0];
finalSignalTab = zeros(lengthInput,1); % Tab of the final processed signal
ecgSignalTab = zeros(lengthInput,1);

while k <= lengthInput
    ecgwindow = ecgSignal(k:min(k+window_size-1, lengthInput)); % Ecg window
    len = length(ecgwindow);
    
    % If window size lower than selected window size
    if len < window_size/2
        k = k-(window_size - len);
        ecgwindow = ecgSignal(k:min(k+window_size-1, lengthInput));
        len = length(ecgwindow);
    end
    
    % Select ecgwindow, ppgwindow, sqi windows and append 0s in the
    % beginning and prepend 0s at the end
    ecgwindow = [zeros(1,Fs) ecgwindow zeros(1,Fs)];
   
    ecg_sig = [0];
    
    if(Fs == 250) % MIT dataset is sampled at 250 Hz not 125
        levelOfDecompostion = 4;
        correctLength = 16; % Signal must be divisible by 2^(Level of decomposition)
    else
        levelOfDecompostion = 3;
        correctLength = 8;
    end

    correctLength = 16; % Signal must be divisible by 2^(Level of decomposition)
    if mod(length(ecgwindow), correctLength)~=0 % Adjusting length for wavelet transform
        ecgwindow = [zeros(1,floor((correctLength-mod(length(ecgwindow),correctLength))/2)) ecgwindow zeros(1,ceil((correctLength-mod(length(ecgwindow),correctLength))/2))];
    end
    
    % ECG Wavelet transform
    [~,swd] = swt(ecgwindow, levelOfDecompostion, 'bior6.8');  % stationary wavelet tranform
    ecg_sig = swd(levelOfDecompostion-1,:)+swd(levelOfDecompostion,:);

    ecg_sig = ecg_sig(round((length(ecg_sig)-len)/2+1):round((length(ecg_sig)+len)/2)); % adjust length to window length
    
% ADDITIVE FUSION (ECG only)
   ecgfin_sig = ecg_sig;
   ecgWavWindow = [ecgWavWindow ecg_sig];
    
    final_sig = ecgfin_sig;
    
    ecgSignalTab(k:k+length(final_sig)-1) = ecgfin_sig;
    finalSignalTab(k:k+length(final_sig)-1) = final_sig; % Tab of signal used for peak detection
       
% Peak detection
    final_sig(1:10) = 0; % This is to ignore the peak in wavelet sequence at start of signal
    [pksmax,locsmax] = findpeaks(final_sig,'MINPEAKDISTANCE',round(0.3*Fs)); % find peaks every 0.3 s
    Vmax = 0.5*mean(pksmax); % set a threshold based on the peaks
    maximatokeep = pksmax>=Vmax; % only peaks above the threshold
    locsmaxima = locsmax([maximatokeep]); % locations of selected peaks
    maxmin_locs=sort(locsmaxima); % sort the locations
    maxmin_val=final_sig(maxmin_locs); % find the peaks at the locations 
    L = 1;
    clear hbloc
    for i = 1:length(maxmin_val)
        if(maxmin_val(i) >= Vmax)
            hbloc(L) = k+maxmin_locs(i)-1; % correct the peak annotation value within this window to that of the signal as a whole
            L = L + 1;
        end
    end
    if L>1
        annotation = [annotation hbloc];
    end
    
    k=k+window_size;
end

% Extract the detected peaks
annotation=annotation(3:end);
annotation = unique(annotation);
% Check for false negatives and positives and correct

peakLocations = [0];

for k = 1:length(annotation)-2
    if peakLocations(end) >= annotation(k)% In case the position at the last position of secondaryLocations is greater than that of annotation
        continue
    end
    rrinterval=(annotation(k+1)-annotation(k))/Fs; % calculate rrinterval
    t=0;
    if rrinterval > 1.2 % if rrinterval is greater than 1.2 s, means missed beat
        signalwindow=finalSignalTab(annotation(k)+round(rrinterval*Fs/3):annotation(k+1)-round(rrinterval*Fs/3));
        if sum(signalwindow)~=0
            [~,idx_val]=max(signalwindow);
            hbloc=annotation(k)+idx_val+round(rrinterval*Fs/3);
            peakLocations=[peakLocations annotation(k) hbloc];
            t=1;
        end  
    end
    if rrinterval<0.3 % if rrinterval is less than 0.4s,means false beat
        if rrinterval<0.1 % if two peaks within 0.1 s, take average
            hbloc=annotation(k+1);
            peakLocations=[peakLocations hbloc];
            t=1;
        else
            rrinterval2=(annotation(k+2)-annotation(k))/Fs; % check interval between the next to next beat
            if rrinterval2 < 1.2 % if less than 1.2, omit k+1 location
                peakLocations=[peakLocations annotation(k) annotation(k+2)];  
                t=1;
            else
                peakLocations=[peakLocations annotation(k)]; % else keep it for next round
                t=1;
            end
        end
    end
    if t==0 % if no corrections made, just add at the end
       peakLocations=[peakLocations annotation(k)]; 
    end
end

peakLocations = [peakLocations annotation(k+1)];
peakLocations = peakLocations(2:end);
ecgWavWindow = ecgWavWindow(2:end);
signal2plot = [ecgSignal]'; % For plotting purpose

if(length(ecgWavWindow) > length(finalSignalTab))
    ecgWavWindow = ecgWavWindow(1:length(finalSignalTab));
end

% Finally get rr interval
RRIntervalSet = zeros(1, length(peakLocations)-1);
for k = 1:length(peakLocations)-1
    rrinterval = (peakLocations(k+1) - peakLocations(k))/Fs; % calculate rrinterval
    RRIntervalSet(k) = rrinterval;
end