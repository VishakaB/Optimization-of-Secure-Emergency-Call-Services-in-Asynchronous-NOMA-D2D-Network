%date: 13 july 2022
%goal: optimize the secrecy 
%for irc
close all
%clc
%clear all

%plotting 
close all
figure (2)
transmit_snrdb_vec = linspace(1,200,20);
for cases = 2:4

C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colros.
%% performance analysis
%secrecy capacity vs SINRd
hold on;
grid on;
plot(transmit_snrdb_vec,smooth(smooth(smooth(cs(cases,:)))),...
    'color',C{cases},'marker','o');
hold on;
plot(transmit_snrdb_vec,smooth(smooth(smooth(wtcs(cases,:)))),...
    'color',C{cases},'LineStyle','--');
xlabel('SNR legitimate');
ylabel('Secrecy capacity');

end

legend({'Proposed:\alpha_{{e}} = 0.01';'Conventional: \alpha_{{e}} = 0.01';...
    'Proposed:\alpha_{{e}} = 0.1';'Conventional: \alpha_{{e}} = 0.1';...
    'Proposed:\alpha_{{e}} = 0.5';'Conventional: \alpha_{{e}} = 0.5';...
    },'Fontsize',9);

