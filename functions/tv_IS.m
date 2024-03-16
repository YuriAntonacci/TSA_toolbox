function ret=tv_IS(A,Su)
%% computation of time varying IS from TV-AR model parameters
% Version for monovariate time series (computes Entropy and information Storage)
% Instead of using the data, uses the VAR parameters (either estimated from data, or theoretical)
% Works under the linear Gaussian assumption
% input:
% A--> AR parmeters (1 x p x N)
% Su--> Residual variance (1 x 1 x N)
% output:
% ret.Sy --> Estimation of the variance of process Y
% ret.Hy --> Entropy of the process Y
% ret.Hy_y --> entropy of Yn given its past
% ret.IS --> time-varying Information storage

%% COVARIANCE MATRICES
M = size(A,1); %number of elements in the system
p=floor(size(A,2)/M); %number of lags in MVAR model
N=size(A,3); % time points
% Obtain F and Delta
Im=eye(M*p);

for n=1:N
    F=[A(:,:,n);Im(1:end-M,:)];% this is A^p Eq.(17)
    Delta=zeros(p*M,p*M); %(this is actually Sigma^p, but use Delta for clarity in code)
    Delta(1:M,1:M)=Su(:,:,n);
    
    % Obtain R_o^p=BigSigma solving the Lyapunov equation: BigSigma = F * BigSigma * F^T + Delta
    BigSigma=dlyap(F,Delta);
    
    %%% Entropy of Yn
    Sy=BigSigma(1,1);
    Hy=real(0.5*log(Sy)+0.5*log(2*pi*exp(1)));
    %%% Conditional entropy of Y given its past
    Sy_Ypast=Su(:,:,n); % variance of innovation is defined as the partial variance
    Hy_Ypast=0.5*log(Sy_Ypast)+0.5*log(2*pi*exp(1));
    
    SY(n,1)=Sy;
    SY_Ypast(n,1)=Sy_Ypast;
    HY(n,1)=Hy;
    HY_Ypast(n,1)=Hy_Ypast;
    % Information storage
    IS(n,1)=Hy-Hy_Ypast;
    
end

ret.Sy=SY;
ret.Sy_y=SY_Ypast;
ret.Hy=HY;
ret.Hy_y=HY_Ypast;
ret.IS=IS;

end