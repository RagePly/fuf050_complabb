
function [a]=getk5_ba_data()
n=1024;
b=importdata('background_97h.txt');
channelnr=linspace(0,n,n);
data=sum(reshape(b(:,2),8192/n,n));
a=[channelnr;data]';
