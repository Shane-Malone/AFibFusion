%[signal,Fs,tm]=rdsamp('3001262_0016.dat',[],[],3); % select record to use for template generation
Fs = 125;
[signal,~,~]=rdsamp('Datasets\p016684\p016684-2188-01-29-00-06.hea',[],3235*Fs,3215*Fs-20,1);
data=double(signal(:,5)); %select the lead used
data=(data-mean(data))/std(data); 
start_seg=data(185:263); %select a window with a single cycle of the signal
start_seg_copy=start_seg;
count=1;
M=39; % length of window= 2M
plot(start_seg)
hold on

for i=632:78:1605 %START: length of window: possible stop point
    segment=data(i-M:i+M,1);
    [acor,lag] = xcorr(segment,start_seg_copy);
    [~,I] = max((acor));
    lagDiff = lag(I);
    circ_segment=circshift(segment,-lagDiff);
    [~,i1]=max(circ_segment);
    [~,i2]=max(start_seg_copy);
    if abs(i1-i2)<2
        plot(circ_segment)
        hold on
        start_seg=start_seg+circ_segment;
        count=count+1;
    end
    if count>10
        break
    end
end
template=start_seg/count;
figure
plot(template)