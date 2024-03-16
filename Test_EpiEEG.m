close all
clear all
clc
%% Analysis of Time-specific information storage during right whisker stimulations
% Performed only on animal IC070523

% load dataset [time steps x EEG channels x Trials]
load('IC070523.mat');
F_f=0.98; %forgetting factor
q=10; % memory of the process 

for epoch=1:size(TOT,3) %epoch cycle
    DATA=TOT(:,:,epoch)'; %Nch x Nsamp
    % TV and local estimation of IS
    for ch=1:size(DATA,1) %channels cycle
        indS=1;
        indE=360;
        % z-score of the data
        DD=zscore(DATA(ch,indS:indE),0,2);
        % model order estimation
        [pottaic,pottmdl] = mos_idMVAR(DD,20,0);
        data_EEG(1,:,:)=DD;
        % identification procedure
        [eAm,eSu]=idMVAR(DD,pottmdl,0);
        % recursive identification of AR model
        [A_e,Su_e]=RLS_ID_AR1_IC(data_EEG,pottmdl,F_f,eAm,eSu);
        % time-varying estimation of IS
        ret_e=tv_IS(A_e,Su_e);
        IS_tv(:,ch,epoch)=ret_e.IS;
        % estimation of local IS 
        lvar = localVAR(eAm,eSu,q);
        out = localInfoStorage(DD',lvar);
        IS_loc(:,ch,epoch) = out.s_y;
        clear out lvar eAm eSu ret_e A_e Su_e data_EEG DD
    end
    clearvars -except IS_tv IS_loc q F_f epoch List_Ch TIME TOT indS indE
end

%% plot of the results over the contralateral and ipsilateral hemisphere
%% Contralateral
% time-varying IS
IS_tv_M=squeeze(nanmean(IS_tv,3));

CH=[10,12,14];

figure
for i=1:length(CH)
    plot(TIME(indS:indE)',IS_tv_M(:,CH(i)));
    hold on
end
legend(List_Ch{CH})
grid on
ylabel('S_{y,n}[nats]')
xlabel('Time [s]')
title('Contralateral')
% local IS
IS_loc_M=squeeze(nanmean(real(IS_loc),3));
figure
for i=1:length(CH)
    plot(TIME(indS:indE)',IS_loc_M(:,CH(i)));
    hold on
end
legend(List_Ch{CH})
grid on
ylabel('s_{y,n}[nats]')
xlabel('Time [s]')
title('Contralateral')
%plot of SEPs
TOT_M=squeeze(nanmean(TOT,3));
figure
for i=1:length(CH)
    plot(TIME',TOT_M(:,CH(i)));
    hold on
end
legend(List_Ch{CH})
grid on
ylabel('Amplitue [V]')
xlabel('Time [s]')
title('Contralateral')
%%  Ipsilateral
% TV-IS
CH=[6,4,2];
figure
for i=1:length(CH)
     plot(TIME(indS:indE)',IS_tv_M(:,CH(i)));
    hold on
end
legend(List_Ch{CH})
grid on
ylabel('S_{y,n}[nats]')
xlabel('Time [s]')
title('Ipsilateral')
% local IS
figure
for i=1:length(CH)
    plot(TIME(50:indE)',IS_loc_M(50:end,CH(i)));
    hold on
end
legend(List_Ch{CH})
legend(List_Ch{CH})
grid on
ylabel('s_{y,n}[nats]')
xlabel('Time [s]')
title('Ipsilateral')
% plot of SEPs
figure
for i=1:length(CH)
    plot(TIME',TOT_M(:,CH(i)));
    hold on
end
legend(List_Ch{CH})
grid on
ylabel('Amplitue [V]')
xlabel('Time [s]')
title('Ipsilateral')


