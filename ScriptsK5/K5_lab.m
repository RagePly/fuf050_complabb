%% K5 lab data analysis

e1 = 482;
e3 = 976; 
e4 = 1048; 

data = data_save;

ch = data(:,1);
co = data(:,2);

rangemin1 = 340;
rangemax1 = 347;
[peak1,fit1,fitrange1]=fitpeak(rangemin1,rangemax1,ch,co);

rangemin2 = 694;
rangemax2 = 705;
[peak2,fit2,fitrange2]=fitpeak(rangemin2,rangemax2,ch,co);

rangemin3 = 748;
rangemax3 = 756;
[peak3,fit3,fitrange3]=fitpeak(rangemin3,rangemax3,ch,co);


p = polyfit([peak1 peak2 peak3],[e1 e3 e4],1);

e = p(1)*ch+p(2);

plot(e, co, 'b')
hold on

plot(p(1)*fitrange1 + p(2), fit1, 'r')
plot(p(1)*fitrange3 + p(2), fit3, 'r')
plot(p(1)*fitrange2 + p(2), fit2, 'r')

xlabel('Energy [keV]')
ylabel('Counts')

hold off





