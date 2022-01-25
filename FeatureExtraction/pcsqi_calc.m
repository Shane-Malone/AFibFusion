function [ppgsqi]=pcsqi_calc(ppgsignal,signal_samplingrate,ppgcurve,curvefreq)
%%%%function to calculate curve matching sqi
%%%ppgsignal - single channel ppgsignal
%%%%signal_samplingrate - sampling frequency of the ppg signal
%%% ppg curve - curve against which correlation of each heartbeat section
%%% is compared.
%%% curvefreq - sampling rate of the curve
if sum(ppgsignal)==0
    ppgsqi=NaN(length(ppgsignal),1);
    ppg_sqi=NaN(length(ppgsignal),1);
else
    Fs=signal_samplingrate;
    if ~isrow(ppgsignal)
        ppgsignal=ppgsignal';
    end
    ppgsqi=zeros(length(ppgsignal),1);
    ppg_sqi=zeros(length(ppgsignal),1);
    ppgcurve = resample(ppgcurve,Fs,curvefreq);
    window_size=length(ppgcurve);
    ppgsignal2=[fliplr(ppgsignal(1:round(window_size/2))) ppgsignal fliplr(ppgsignal((length(ppgsignal)-round(window_size/2)):length(ppgsignal)))];%%%mirroring signal towards the beginning and end

    k=1;
    x=ppgcurve';
    y=fliplr(x);
    y=[y zeros(1,window_size-1)];
    P=gallery('circul',y);
    P=P(1:window_size,:);
    lag=[-(window_size-1):1:window_size-1];
    window=[ppgsignal2(max(1,k-round(window_size/2)+1):k-1) ppgsignal2(k) ppgsignal2(k+1:min(k+round(window_size/2)-1,length(ppgsignal2)))]; %%%mirroring signal towards the beginning and end
    if k-round(window_size/2)<1
            window=[zeros(1,round(window_size/2)-k) window];
    end
    if k+round(window_size/2)>=length(ppgsignal2)
        window=[window zeros(1,(window_size-length(window)))];
    end
    window1=diag(window);
    tot=window1*P;
    total=sum(tot);
    total1=fliplr(total);
    [~,I]=max(total1);
    beg_last=window(1);
    correction=lag(I);
    corrected_curve=circshift(ppgcurve',-correction);
    diff=var(corrected_curve-window);
    ppg_sqi(k)=diff;
    if ppg_sqi(k)~=0
        ppgsqi(k)=1/ppg_sqi(k);
    else
        ppgsqi(k)=0.014;
    end
    if ppgsqi(k)<0 %%%limiting negative SQI values to 0.0001
        ppgsqi(k)=0;
    end
    


    for k=2:length(ppgsignal2)
        window=window(2:end);
        if round(window_size/2)+k-1>=length(ppgsignal2)
            window=[window zeros(1,(window_size-length(window)))];
        else
            window(window_size)=ppgsignal2(round(window_size/2)+k-1);
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
    	corrected_curve=circshift(ppgcurve',-correction);
        diff=var(corrected_curve-window);
        ppg_sqi(k)=diff;
        ppgsqi(k)=1/ppg_sqi(k);
        ppgsqi(k)=5.9904*log(ppgsqi(k))+53.986; %%%from curve fitting
    
        if ppgsqi(k)<0 %%%limiting negative SQI values to 0.0001
            ppgsqi(k)=0;
        end
    end
    ppgsqi=ppgsqi(round((length(ppgsqi)-length(ppgsignal))/2+1):round((length(ppgsqi)+length(ppgsignal))/2));
end

