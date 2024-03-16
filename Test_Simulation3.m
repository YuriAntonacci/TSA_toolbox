%% Test Simulation3
% Time-specific information storage with innovation variance and AR coefficient 
% which varyi as a periodic square waveform (scenario 3 - Fig.1c main document)
close all;clear all;clc
%% generation of time variant innovation variance and AR coefficient
rng(1,'twister')
cmax      =  0.3;    % max amplitude
fs        =  1000;    % sampling frequency
f_osc     =  0.05;    % frequency of oscillation 
p         =    1;    % model order        

%% Construction of non-stationary VAR time series

nobs = 60000;                        % 10 s of observation (data samples)
k = 1:nobs;                         % vector of time steps
t = (k)/fs;                         % vector of times
DC=50;                              % duty-cycle
c = -cmax*square(2*pi*f_osc*t,DC);  % causal coefficient varying as square periodic waveform
Su_T(1,1,:)=c+cmax*2;               % theoretical vector of innovation variance 
A_T(1,1,:) = c+cmax*2;              % theoretical AR coefficient               
% Data generation
rng(1) % always get the same seed to fix the realization
Y = var_nonstat(A_T,Su_T,1);
% theoretical IS value
ret_t=tv_IS(A_T,Su_T);

%% Time-var IS estimation
fFactor=[0.98];
LIM=[-1, 5];
% RLS algorithm for TV-IS
[A_e,Su_e]=RLS_ID_AR1(Y,p,fFactor);
% comptation of IS
ret_e=tv_IS(A_e,Su_e);

%% Plot of the results
pts=200;
figure
set(gcf, 'Color', 'w');
% derivate of square function
D=abs(diff(c));
[ind]=find(D==max(D));
for pp=1:numel(ind)
    xline(t(ind(pp)),'--k')
end
hold on
plot(t(pts:end),ret_e.IS(pts:end,1),'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
hold on
plot(t(pts:end),[ret_t.IS(pts:end,1)],'k','LineWidth',1);
ylabel('S_{y,n} [nats]')
xlabel('Time [s]')
xlim([0.2 60])
xline(60,'--k')
ylim(LIM);
% plot of average value and variance across time
IND=zeros(1,length(ind)+1);
IND(1)=1;
IND(2:end)=ind;
MEAN=[4];
VAR=[-0.5];
for oo=1:length(IND)-1
    M=nanmean(ret_e.IS(IND(oo):IND(oo+1)));
    V=nanvar(ret_e.IS(IND(oo):IND(oo+1)));
    text(t(round(IND(oo+1)-((IND(oo+1)-IND(oo))/2))),MEAN,sprintf('%4.3f',round(M,3,'Significant')) ,'FontSize',...
        10,'Color','k','FontName','TimesNewRoman','HorizontalAlignment','center');
    text(t(round(IND(oo+1)-((IND(oo+1)-IND(oo))/2))),VAR,sprintf('(%4.3f)',round(V,3,'Significant')),'FontSize',...
        10,'Color','k','FontName','TimesNewRoman','HorizontalAlignment','center');
end
set(gcf,'units','centimeters','position',[0,0,12,12])



%% local information storage estimation and plot of the results
Y=squeeze(Y)';
[eAm,eSu]=idMVAR(Y,p,0);
lvar = localVAR(eAm,eSu,10);
out = localInfoStorage(Y',lvar);
is_total = out.s_y;

LIM=[-15, 15];
figure
set(gcf, 'Color', 'w');
xlim([0.2 60])
ylim(LIM)

MEAN=[12];
VAR=[-12];

for oo=1:length(IND)-1
    M=nanmean(is_total(IND(oo):IND(oo+1)));
    V=nanvar(is_total(IND(oo):IND(oo+1)));
    text(t(round(IND(oo+1)-((IND(oo+1)-IND(oo))/2))),MEAN,sprintf('%4.3f',round(M,3,'Significant')) ,'FontSize',...
        10,'Color','k','FontName','TimesNewRoman','HorizontalAlignment','center');
    text(t(round(IND(oo+1)-((IND(oo+1)-IND(oo))/2))),VAR,sprintf('(%4.3f)',round(V,3,'Significant')),'FontSize',...
        10,'Color','k','FontName','TimesNewRoman','HorizontalAlignment','center');
end
hold on
for pp=1:numel(ind)
    xline(t(ind(pp)),'--k')
end
hold on
plot(t(pts:end),is_total(pts:end,1),'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
hold on
plot(t(pts:end),[ret_t.IS(pts:end,1)],'k','LineWidth',1);
ylabel('s_{y,n} [nats]')
xlabel('Time [s]')
set(gcf,'units','centimeters','position',[0,0,12,12])


