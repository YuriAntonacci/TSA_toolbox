%% computation of local Information Storage from covariance matrices

function out = localInfoStorage(Y,ret)

%   Computes Local Information Dynamics analytically for a stationary mvar(p) process:
%
%   INPUT: 
%   time series Y (one time series in column)
%   ret structure with
%   S_Yn  -   variance of the process
%   S_Ynq  -   covariance matrix of the past state of the process (Wn)
%   S_YnYnq  - joint covariance matrix between the present and past states of the process [Yn, Wn]
%   OUTPUT:
%   out structure with local Information Storage

S_Yn=ret.S_Yn;
S_Ynq=ret.S_Ynq;
S_YnYnq=ret.S_YnYnq;

q=size(S_Ynq,1);
N=size(Y,1);

iS_Ynq=inv(S_Ynq);
iS_YnYnq=inv(S_YnYnq);

s_y_p=nan*ones(N,1); s_y_n=s_y_p;
for n=q+1:N
    yn=Y(n);    % present state
    ynq=Y(n-1:-1:n-q)'; % past state
    s_y_p(n)=0.5*(yn^2)/S_Yn + 0.5*ynq*iS_Ynq*ynq';
    s_y_n(n)=-0.5*[yn ynq]*iS_YnYnq*[yn ynq]';       
end
s_y=s_y_p+s_y_n;

out.s_y=s_y + ret.Sy;
out.Sy = ret.Sy;

out.s_y_p=s_y_p;
out.s_y_n=s_y_n;
