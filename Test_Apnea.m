clear;clc; close all
%% Analysis of Time-specific information storage during sleep apnea

load('apnea.mat'); % x=RR, y=RESP
RESP=zscore(y,0,1);
ats=[43 121 202 289 380 488 547 618 681 778 875 972 1040 1119]; %apnea tstart
ate=[82 168 262 355 436 521 579 655 752 841 928 1006 1089 1172];  %apnea tend

ats1=[1  43 83  121 169 202 263 289 356 380 437 488 522 547 580 618 656 681 753 778 842 875 929 972  1007   1040 1090 1119 1173]; 
ate1=[42 82 120 168 201 262 288 355 379 436 487 521 546 579 617 655 680 752 777 841 874 928 971 1006 1039 1089 1118 1172 1200];  

numsurr = 100;
lo_loc = 0; lh_loc = 100;
lo_tv = 0; lh_tv = 100;
pmax = 20;

for is = 1:numsurr+1 %cycle on surrogates
    % random shuffling procedure
    if is ~=1
        RESP = surrshuf(RESP);
    end
    % model order estimation
    [pottaic,pottmdl,aic,mdl] = mos_idMVAR(RESP',pmax,0);
    % identification procedure
    [eAm_R,eSu_R]=idMVAR(RESP',pottmdl,0);
    
    % local- IS estimation
    q=10;
    lvarx = localVAR(eAm_R,eSu_R,q);
    out_R = localInfoStorage(RESP,lvarx);
    ISl_R(is,:) = out_R.s_y;

    % RLS - Initial Conditions
    R=permute(RESP,[3,2,1]);
    [A_x,Su_x]=RLS_ID_AR1_IC(R,pottmdl,0.98,eAm_R,eSu_R);
    ret_x=tv_IS(A_x,Su_x);
    IStv_R(is,:)=ret_x.IS;
    
end
ISl_th=prctile(ISl_R(2:end,:),[lo_loc lh_loc]);
IStv_th=prctile(IStv_R(2:end,:),[lo_tv lh_tv]);

%% Figure
fsax=18; fstx=20;
FontSize = 15;
set(0,'defaultAxesFontSize', fsax);
set(0,'defaultTextFontSize', fstx);
set(0,'defaultTextFontName', 'Times');
set(0,'defaultAxesFontName', 'Times');

lw=1.1; yminy=-4; ymaxy=6; yminx=-2.5; ymaxx=3;
col_as=[0.5 0.5 0.5];
col_ae=[0.8 0.2 0.2];
col1 = [0,139,139]/255;
col2 = [255,99,71]/255;
colorsh=[0.88 0.88 0.88];

t=(1:1:length(RESP))';
fh=figure(50);
set(gcf, 'Position', [50 50 1000 600]);
set(gcf, 'Color', 'w');
sfh1=subplot(3,1,1,'Parent',fh);
for i=1:length(ats)
    hu=fill([ats(i) ats(i) ate(i) ate(i)]',[yminy ymaxy ymaxy yminy]',colorsh); hold on
    set(hu,'edgecolor','white');
end
plot(y,'color','k','linewidth',lw); ylim([yminy ymaxy]);
ylabel('RESP')
set(gca, 'YTick', [])
sfh1.Position = sfh1.Position + [0 0.02 0 0.037];

% local IS
sk=2*q;
t=(sk:1:length(y));
ymin=-11; ymax=10;
sfh3=subplot(3,1,2,'Parent',fh);
for i=1:length(ats)
    hu=fill([ats(i) ats(i) ate(i) ate(i)]',[ymin ymax ymax ymin]',colorsh); hold on
    set(hu,'edgecolor','white');
end
p1 = patch([q+1:1200,1200:-1:q+1]',[ISl_th(1,q+1:end),fliplr(ISl_th(2,q+1:end))]',col1);
set(p1,'facecolor',col1,'edgecolor','none','facealpha',0.5);
plot(1:1200,ISl_R(1,:),'color',col1,'linewidth',lw);

for i=1:length(ats1)
    M=nanmean(ISl_R(1,ats1(i):ate1(i)));
    V=nanvar(ISl_R(1,ats1(i):ate1(i)));
    MM_loc(i,1)=M;
    VV_loc(i,1)=V;
end

ylim([ymin ymax])
ylabel('s_{RESP,n} [nats]')
sfh3.Position = sfh3.Position + [0 -0.03 0 0.037];

% time-varying IS
ymin=-0.1; ymax=2.6;
sfh4=subplot(3,1,3,'Parent',fh);
for i=1:length(ats)
    hu=fill([ats(i) ats(i) ate(i) ate(i)]',[ymin ymax ymax ymin]',colorsh); hold on
    set(hu,'edgecolor','white');
end
p2 = fill([q+1:1200,1200:-1:q+1]',[IStv_th(1,q+1:end), fliplr(IStv_th(2,q+1:end))]',col2);
set(p2,'edgecolor','none','facealpha',0.5);
plot(1:1200,IStv_R(1,:),'color',col2,'linewidth',lw);

for i=1:length(ats1)
    M=nanmean(IStv_R(1,ats1(i):ate1(i)));
    V=nanvar(IStv_R(1,ats1(i):ate1(i)));
    MM_tv(i,1)=M;
    VV_tv(i,1)=V;
end

ylim([ymin ymax])
ylabel(' S_{RESP,n} [nats]')
sfh4.Position = sfh4.Position + [0 -0.03 0 0.037];
fh.WindowState = 'maximized'; 