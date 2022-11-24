%date: 13 july 2022
%goal: optimize the secrecy 
%for irc
close all
clc
clear all

%% initialization 
K = 3;
max_dist     = 100;%meters
max_eta      = 10;
etath        = 4;%change this 
noisepower   = 0.1;
max_tx_power = 500;%change this
B            = 1;%channel bandwidth
timeslot     = 1;
N            = 10^3;

transmit_snrdb_vec = linspace(1,100,20);%practical range is 0-100db
total = length(transmit_snrdb_vec);
t = 1;%sinr index 
cs = zeros(K,total);
wtcs = zeros(K,total);
alpha1_vec = [0.0001,0.001,0.1,0.5];
alpha2_vec = [0.4,0.3,0.15,0.1];
receive_pow_ratio_vec =linspace(0,10,4);%change here

for cases = 1:4
t = 1; 
receive_pow_ratio = receive_pow_ratio_vec(2);%?

for transmit_snrdb = transmit_snrdb_vec
  
time_offset = 0.15;

%Distances of users from rx
dist_k = max_dist*abs(randn(K,1));

%sorted transmit power vector %descending 
dist_vec = sort(dist_k,'ascend'); 

% Path loss exponent
eta_k   = max_eta*abs(randn(K,1));
eta_vec = sort(eta_k,'ascend'); 

% unsorted transmit power vector
transmitpow_k = max_tx_power*abs(randn(K,1));

%sorted transmit power vector %descending 
power_vec = sort(transmitpow_k,'descend'); 
power_vec_e = sort(transmitpow_k,'descend'); 
power_vec(1) = 10^(transmit_snrdb/10)*noisepower;
%receive_pow_ratio = receive_pow_ratio_vec(receive_pow_ratioi);
%change here
for d = 2: K
    power_vec(d) = power_vec(d-1)/10^(receive_pow_ratio);
end

power_vec_e(1) = (alpha1_vec(cases))*10^(transmit_snrdb/10)*noisepower;

for d = 2: K
    power_vec_e(d) = power_vec_e(d-1)/10^(receive_pow_ratio);
end
%tx power vec
%power_vec = tx_power_percentage_vec.*tx_pow_k;
eta_vec = 1.5;
pathloss_exp = sqrt(dist_vec.^-eta_vec);

%channel coefficients of each user vec
h_vec =  pathloss_exp.*sqrt(power_vec/2).*(randn(1,N)+1i*randn(1,N))/sqrt(2);
h_vec_e =  pathloss_exp.*sqrt(power_vec_e/2).*(randn(1,N)+1i*randn(1,N))/sqrt(2);

%channel gains of each user vec
g_vec = (abs(h_vec)).^2;
g_vec_e = (abs(h_vec_e)).^2;

for k = 1:K
    nsym(k,1) = K-(k-1);
    user_strength(k,1) = k;
end

delta_mat = zeros(K,K);
delta_mat(1,:) = zeros(K,1);
delta_mat(2:K,:) = time_offset*ones(K-1,K);%B1,... Bn, C1....,Cn, ....... %Z1,....Zn

reverse_delta_mat(K,:)  = zeros(K,1);
reverse_delta_mat(1:K-1,:) = time_offset*ones(K-1,K);%A1, A2

desired_id = 1;%initiation
for j = 1:K%interference vector loop
for k = 1:K%interfering user
    %interference vec %only from the next neighbor user
    if k ~= desired_id & k == desired_id+1 & desired_id == 1        
        sumsym_dur_vec(desired_id,1) = sum(delta_mat(desired_id+1,:))       
    elseif k ~= desired_id & k == desired_id+1 & desired_id > 1 & desired_id <K
        sumsym_dur_vec(desired_id,1) = sum(delta_mat(desired_id+1,:))+...
            sum(reverse_delta_mat(desired_id-1,:))
    elseif desired_id==K
        sumsym_dur_vec(desired_id,1) = sum(reverse_delta_mat(desired_id-1,:));       
    end
 desired_id = desired_id+1;
end
end
sumsym_dur_vec

interf_vec = zeros(K,1);
desired_id = 1;
for j = 1:K%interference vector loop
for k = 1:K%neighbor users index
     if k ~= desired_id  %strongest
        interf_vec(desired_id) = interf_vec(desired_id)+sum(power_vec(k).*...
            mean(g_vec(k,:),2).*delta_mat(desired_id,:));
     elseif desired_id == 1
        if(k<K)
        interf_vec(desired_id) = sum(power_vec(k+1).*...
            mean(g_vec(k+1,:),2).*delta_mat(desired_id+1,:));
        end
     end
end
desired_id = desired_id+1;
end

%% optimization
%capacity of legitimate user 
sinr_thmin = 1e-3;
sinr_thmax = 1;
ratemin = 1e6;
ratemax = 1e9;
B = 1e6;

SINRd = sum(power_vec(1:K).*mean(g_vec(1:K,:),2)./...
    (noisepower^2 + (interf_vec)));
SINRe = sum(power_vec_e(1:K).*mean(g_vec_e(1:K,:),2)./...
    (noisepower^2 + (interf_vec)));

for k = 1:K
    K_vec(k,1) = K-(k-1);
end
for k = 1:K
    e_vec(k,1) = K-(k-1);
end

cvx_begin quiet
   variable decision_uk(K,1)
   dual variables var1 var2 var3 var4
   minimize(-decision_uk'*K_vec)
   subject to
         var1 : -sum(decision_uk)+ sum(decision_uk.^2)<=0;
         var2 : decision_uk.*((noisepower^2 + interf_vec +...
          power_vec(1:K).*mean(g_vec(1:K,:),2))...
                -power_vec(1:K).*mean(g_vec(1:K,:),2)*...
                (1+1/sinr_thmin))- (decision_uk.*((noisepower^2 + interf_vec +...
          power_vec_e(1:K).*mean(g_vec_e(1:K,:),2))...
                -power_vec_e(1:K).*mean(g_vec_e(1:K,:),2)*...
                (1+1/sinr_thmax))) <= 0
         var3: decision_uk >=0 
         var4: decision_uk(1:K-1)>= decision_uk(2:K) 
cvx_end

%decision_uk = decision_uk>0.25;

cs(cases,t) = (log((1+sum(power_vec(1:K).*mean(g_vec(1:K,:),2)./...
    (noisepower^2 + (interf_vec)))))...
       -log(1+sum(decision_uk.*power_vec_e(1:K).*mean(g_vec_e(1:K,:),2)./...
       (noisepower^2 + (interf_vec)))));%maximize secrecy
wtcs(cases,t) = (log((1+sum(power_vec(1:K).*mean(g_vec(1:K,:),2)./...
    (noisepower^2 + (interf_vec)))))...
       -log(1+sum(power_vec_e(1:K).*mean(g_vec_e(1:K,:),2)./...
       (noisepower^2 + (interf_vec)))));
t = t+1;
end
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

legend({'opt \alpha_{{e}_1} = 0.001';'wt opt \alpha_{{e}_1} = 0.001';...
    'opt \alpha_{{e}_2} = 0.01';'wt opt \alpha_{{e}_2} = 0.01';...
    'opt \alpha_{{e}_3} = 0.1';'wt opt \alpha_{{e}_3} = 0.1';...
    'opt \alpha_{{e}_4} = 0.5';'wt opt \alpha_{{e}_4} = 0.5'});