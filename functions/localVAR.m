%% computation of global information storage and covariance matrices from AR model parameters
%% univariate version of localCOV.m

function ret = localVAR(Am,Su,q)

%   Computes Information Dynamics analytically for a stationary var(p) process:
%   X_n=A_1*X_{n-1}+A_2*X_{t-n}+...+A_p*X_{n-p}+E_n
%
%   INPUT: 
%   Am  -   generalized connectivity matrix A=(A_1 A_2 ... A_p)
%   Su  -   covariance matrix for E_n
%   q   -   number of lags
%   OUTPUT:
%   ret structure with covariance matrices and global Information Dynamics

p=size(Am,2); %number of lags in MVAR model

R=NaN*ones(1,1,q+1); % prepare covariance matrices, (:,:,1) is lag 0, (:,:,q+1) is lag q 

% Obtain F and Delta
Im=eye(p);
F=[Am;Im(1:end-1,:)];% this is A^p
Delta=zeros(p,p); %(this is actually Sigma^p, but use Delta for clarity in code)
Delta(1,1)=Su;

% Obtain R_o^p=BigSigma solving the Lyapunov equation: BigSigma = F * BigSigma * F^T + Delta
BigSigma=dlyap(F,Delta);

% extract correlation for lag 0... p-1: R(0),...,R(p-1)
for i=1:p
    R(:,:,i)=BigSigma(1,i);
end

% Yule-Walker solution  for lags >= p
for k=p+1:q+1
    Rk=R(:,:,k-1:-1:k-p);
    Rm=[];
    for ki=1:p
        Rm=[Rm; Rk(:,:,ki)];
    end
    R(:,:,k)=Am*Rm;
end
Ry=squeeze(R)';

S_Yn=Ry(1);     % innovation variance
S_Ynq=toeplitz(Ry(1:q),Ry(1:q));    % covariance matrix of Wn 
S_Yn_Ynq=Ry(2:q+1);         % Cross-covariance matrix between Yn and Wn
% S_Ynq_Yn=S_Yn_Ynq';
S_YnYnq=[[S_Yn S_Yn_Ynq];[S_Yn_Ynq' S_Ynq]];    % joint covariance matrix between present Yn, and past state Wn of the process

Sy=0.5*log((S_Yn*det(S_Ynq))/det(S_YnYnq));     % Information storage

ret.Sy=Sy;
ret.S_Yn=S_Yn;
ret.S_Ynq=S_Ynq;
ret.S_YnYnq=S_YnYnq;
ret.BigSigma=BigSigma;

