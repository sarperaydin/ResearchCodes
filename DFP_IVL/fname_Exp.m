clear
clc
a=2;
N=5;


FileName=['Agents',num2str(N),'-',datestr(now,'mmmmddyyyy-HHMMAM'),'.mat'];

%save 'proj'+datestr(now, 'dd-mmm-yyyy') a;

save(FileName,'a');