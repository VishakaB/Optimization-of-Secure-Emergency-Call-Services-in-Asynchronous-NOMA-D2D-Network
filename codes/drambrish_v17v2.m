%date: 13 july 2022
%goal: optimize the secrecy 
%for irc
close all
%clc
%clear all

%plotting 

for cases = 2:4

C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colros.
%% performance analysis
%secrecy capacity vs SINRd
hold on;
grid on;
plot(transmit_snrdb_vec,smooth(smooth(smooth(cs(cases,:)))),'color',C{cases},'marker','o');
hold on;
plot(transmit_snrdb_vec,smooth(smooth(smooth(wtcs(cases,:)))),'color',C{cases},'marker','*');
xlabel('SNR legitimate');
ylabel('Secrecy capacity');

end
proposed = smooth(smooth(smooth(cs(cases,:))));
conv = smooth(smooth(smooth(wtcs(cases,:))));
save proposed.mat;
save conv.mat;

legend({...
%{'opt \alpha_{{e}_1} = 0.001';'wt opt \alpha_{{e}_1} = 0.001';...
    'proposed: \alpha_{{e}} = 0.01';'conventional: \alpha_{{e}} = 0.01';...
    'proposed \alpha_{{e}} = 0.1';'conventional: \alpha_{{e}} = 0.1';...
    'proposed: \alpha_{{e}} = 0.5';'conventional: \alpha_{{e}} = 0.5'},'Fontsize',9);
 