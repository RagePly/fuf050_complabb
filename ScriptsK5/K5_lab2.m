%% Kurieplot

data = data_Barium_save;
ba_data = getk5_ba_data();

%data(:,2) = data(:,2) - ba_data(:,2);
e = p(1)*data(:,1)+p(2);
co = data(:, 2);

rangemin1 = 613;
rangemax1 = 634;
[peak1,fit1,fitrange1]=fitpeak(rangemin1,rangemax1,e,co);

rangemin2 = 649;
rangemax2 = 663;
[peak2,fit2,fitrange2]=fitpeak(rangemin2,rangemax2,e,co);

figure(101)
semilogy(e, data(:,2))

[kuriedata] = kurieplot(e, data, 56);

figure(102)
semilogy(e,kuriedata)

[slope,offset]=kuriefit(e,kuriedata,700:1000); 

Q1 = -offset/slope; %high energy

[slope2,offset2]=kuriefit(e,kuriedata,270:450);

Q2 = -offset2/slope2; % low energy

f1=log10f(56,Q1);
f2=log10f(56,Q2);



