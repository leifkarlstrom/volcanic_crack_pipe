clear all
close all

x=0:1:1000;

%BG = 'henrys_law';
BG = {'specified_n','henrys_law','Halemaumau_mean','Halemaumau_max'};
%BG = 'Halemaumau_mean';

for ii=1:length(BG)

[rho K c a b, p, n]=magma_st(x,BG{ii});

subplot(1,3,1)
plot(n,x)
xlabel('n')
hold on

subplot(1,3,2)
semilogx(rho,x)
hold on
xlabel('rho')
xlim([10,4000])

subplot(1,3,3)
plot(c,x)
hold on
xlabel('c')
end
legend(BG)




