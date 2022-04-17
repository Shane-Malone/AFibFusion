% function to calculate curve matching sqi for ECG
% Ecgsignal - single channel Ecgsignal
% csqi_calc returns the csqi value
% Author - Arlene John, Shane Malone
% Initial code developed by Arlene. Adapted for use in fusion algorithm by
% Shane Malone.

% signal_samplingrate - sampling frequency of the ecg signal
% ecg curve - curve against which correlation of each heartbeat section
% is compared.
% curvefreq - sampling rate of the curve

function [ecgsqi]=csqi_calc(Ecgsignal,signal_samplingrate,ecgcurve,curvefreq)

if sum(Ecgsignal)==0
    ecgsqi=NaN(length(Ecgsignal),1);
    ecg_sqi=NaN(length(Ecgsignal),1);
else
    Fs=signal_samplingrate;
    if ~isrow(Ecgsignal)
        Ecgsignal=Ecgsignal';
    end
    ecgsqi=zeros(length(Ecgsignal),1);
    ecg_sqi=zeros(length(Ecgsignal),1);
    ecgcurve = resample(ecgcurve,Fs,curvefreq);
    window_size=length(ecgcurve);
    Ecgsignal2=[fliplr(Ecgsignal(1:round(window_size/2))) Ecgsignal fliplr(Ecgsignal((length(Ecgsignal)-round(window_size/2)):length(Ecgsignal)))];% mirroring signal towards the beginning and end

    k=1;
    x=ecgcurve';
    y=fliplr(x);
    y=[y zeros(1,window_size-1)];
    P=gallery('circul',y); % Generate circulant matrix(Toeplitz)
    P=P(1:window_size,:);
    lag=[-(window_size-1):1:window_size-1];
    window=[Ecgsignal2(max(1,k-round(window_size/2)+1):k-1) Ecgsignal2(k) Ecgsignal2(k+1:min(k+round(window_size/2)-1,length(Ecgsignal2)))]; % mirroring signal towards the beginning and end
    if k-round(window_size/2)<1
            window=[zeros(1,round(window_size/2)-k) window];
    end
    if k+round(window_size/2)>=length(Ecgsignal2)
        window=[window zeros(1,(window_size-length(window)))];
    end
    window1=diag(window);
    tot=window1*P;
    total=sum(tot);
    total1=fliplr(total);
    [~,I]=max(total1);
    beg_last=window(1);
    correction=lag(I);
    corrected_curve=circshift(ecgcurve',-correction);
    diff=var(corrected_curve-window);
    ecg_sqi(k)=diff;
    if ecg_sqi(k)~=0
        ecgsqi(k)=1/ecg_sqi(k);
    else
        ecgsqi(k)=0.014;
    end
    if ecgsqi(k)<0 % limiting negative SQI values to 0.0001
        ecgsqi(k)=0;
    end
    


    for k=2:length(Ecgsignal2)
        window=window(2:end);
        if round(window_size/2)+k-1>=length(Ecgsignal2)
            window=[window zeros(1,(window_size-length(window)))];
        else
            window(window_size)=Ecgsignal2(round(window_size/2)+k-1);
        end
        last_term=window(end);
        new_row=last_term.*P(end,:);
        first_row=beg_last.*P(1,:);
        total=total-first_row;
        total=circshift(total,-1);
        total=total+new_row;
        beg_last=window(1);
        total1=fliplr(total);
        [~,I]=max(total1);
        correction=lag(I);
    	corrected_curve=circshift(ecgcurve',-correction);
        diff=var(corrected_curve-window);
        ecg_sqi(k)=diff;
        ecgsqi(k)=1/ecg_sqi(k);


        ecgsqi(k)=(2.086*10^-8)*exp(12.35*ecgsqi(k));
    

        if ecgsqi(k) < 0
             ecgsqi(k) = 0.0001;
        elseif ecgsqi(k) > 20
            ecgsqi(k) = 20;
        end
    end
    ecgsqi=ecgsqi(round((length(ecgsqi)-length(Ecgsignal))/2+1):round((length(ecgsqi)+length(Ecgsignal))/2));
end


