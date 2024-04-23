close all
clear all
clc
%% Analysis of Time-specific information storage during right whisker stimulations
% Performed only on animal IC070523

% load dataset [time steps x EEG channels x Trials]
load('IC070523.mat');
F_f=0.98; %forgetting factor
indS=1;
indE=360;
for ch=1:size(TOT,2) %channels cycle
    disp(ch)
    % z-score of the data
    DD=zscore(TOT(:,ch),0,1);
    % model order estimation
    [pottaic,pottmdl] = mos_idMVAR(DD',30,0);
    data_EEG(1,1,:)=DD;
    % identification procedure
    [eAm,eSu]=idMVAR(DD',pottmdl,0);
    % recursive identification of AR model
    [A_e,Su_e]=RLS_ID_AR1_IC(data_EEG,pottmdl,F_f,eAm,eSu);
    % time-varying estimation of IS
    ret_e=tv_IS(A_e,Su_e);
    IS_tv(:,ch)=ret_e.IS;
    % estimation of local IS
    lvar = localVAR(eAm,eSu,pottmdl);
    out = localInfoStorage(DD,lvar);
    IS_loc(:,ch) = out.s_y;
    clear out lvar eAm eSu ret_e A_e Su_e data_EEG DD
end
% Reshape data to compute average across trials
length_trial=indE;
for ich = 1:size(IS_tv,2)
    IS_tv_trials(:,:,ich)=reshape(squeeze(IS_tv(:,ich)),[length_trial, size(IS_loc,1)/length_trial]);
    IS_loc_trials(:,:,ich)=reshape(squeeze(IS_loc(:,ich)),[length_trial, size(IS_loc,1)/length_trial]);
    TOT_trials(:,:,ich)=reshape(squeeze(TOT(:,ich)),[length_trial, size(IS_loc,1)/length_trial]);
end
IS_tv=IS_tv_trials;
IS_loc=IS_loc_trials;
TOT=TOT_trials;
%% plot of the results over the contralateral and ipsilateral hemisphere
%% Contralateral
% time-varying IS
IS_tv_M=squeeze(nanmean(IS_tv,2));

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
IS_loc_M=squeeze(nanmean(real(IS_loc),2));
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
TOT_M=squeeze(nanmean(TOT,2));
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


