function [abpsqi]=Acsqi_calc(Abpsignal,signal_samplingrate,abpcurve,curvefreq)
%%%%function to calculate curve matching sqi
%%%Abpsignal - single channel Abpsignal
%%%%signal_samplingrate - sampling frequency of the ecg signal
%%% abp curve - curve against which correlation of each heartbeat section
%%% is compared.
%%% curvefreq - sampling rate of the curve


if sum(Abpsignal)==0
    abpsqi=NaN(length(Abpsignal),1);
    abp_sqi=NaN(length(Abpsignal),1);
else
    Fs=signal_samplingrate;
    if ~isrow(Abpsignal)
        Abpsignal=Abpsignal';
    end
    abpsqi=zeros(length(Abpsignal),1);
    abp_sqi=zeros(length(Abpsignal),1);
    abpcurve = resample(abpcurve,Fs,curvefreq);
    window_size=length(abpcurve);
    Abpsignal2=[fliplr(Abpsignal(1:round(window_size/2))) Abpsignal fliplr(Abpsignal((length(Abpsignal)-round(window_size/2)):length(Abpsignal)))];%%%mirroring signal towards the beginning and end

    k=1;
    x=abpcurve';
    y=fliplr(x);
    y=[y zeros(1,window_size-1)];
    P=gallery('circul',y); % Generate circulant matrix(Toeplitz)
    P=P(1:window_size,:);
    lag=[-(window_size-1):1:window_size-1];
    window=[Abpsignal2(max(1,k-round(window_size/2)+1):k-1) Abpsignal2(k) Abpsignal2(k+1:min(k+round(window_size/2)-1,length(Abpsignal2)))]; %%%mirroring signal towards the beginning and end
    if k-round(window_size/2)<1
            window=[zeros(1,round(window_size/2)-k) window];
    end
    if k+round(window_size/2)>=length(Abpsignal2)
        window=[window zeros(1,(window_size-length(window)))];
    end
    window1=diag(window);
    tot=window1*P;
    total=sum(tot);
    total1=fliplr(total);
    [~,I]=max(total1);
    beg_last=window(1);
    correction=lag(I);
    corrected_curve=circshift(abpcurve',-correction);
    diff=var(corrected_curve-window);
    abp_sqi(k)=diff;
    if abp_sqi(k)~=0
        abpsqi(k)=1/abp_sqi(k);
    else
        abpsqi(k)=0.014;
    end
    if abpsqi(k)<0 %%%limiting negative SQI values to 0.0001
        abpsqi(k)=0;
    end
    


    for k=2:length(Abpsignal2)
        window=window(2:end);
        if round(window_size/2)+k-1>=length(Abpsignal2)
            window=[window zeros(1,(window_size-length(window)))];
        else
            window(window_size)=Abpsignal2(round(window_size/2)+k-1);
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
    	corrected_curve=circshift(abpcurve',-correction);
        diff=var(corrected_curve-window);
        abp_sqi(k)=diff;
        abpsqi(k)=1/abp_sqi(k);
        abpsqi(k)=5.026*log(abpsqi(k))+47.78; %%%from curve fitting
    
        if abpsqi(k)<0 %%%limiting negative SQI values to 0.0001
            abpsqi(k)=0;
        end
    end
    abpsqi=abpsqi(round((length(abpsqi)-length(Abpsignal))/2+1):round((length(abpsqi)+length(Abpsignal))/2));
end


