%% var_nonstat
%
% Generate random multi-trial non-stationary Gaussian VAR time series
% input:
%     A    -->      time-varying VAR coefficients matrix [1 x p x N] (N - total number of time points)
%     Su   -->      time-varying residuals variance [1 x 1 x N]
%     Ntr  -->      number of realizations (default: 1)
%     
% output
%     Y -->         multi-trial Gaussian VAR time series [R x 1 x N]
%     E -->         residuals time series [R x 1 x N]
%
%% Description
%
% Return |R| independent random Gaussian non-stationary VAR time series, with
% time-varying coefficients |A| and residuals covariances |Su|. The last index
% in |A| and |Su| is the time index (the number of observations is taken as the
% size of the last dimension). If |Su| is 2 dimensional, then it is replicated
% at each time.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
%%

function Y = var_nonstat(A,Su,R)

if nargin < 3 || isempty(R), R = 1; end % single trial

[n,p,m] = size(A);
[n2,n3,m1] = size(Su);
assert(n2 == n && n3 == n,'residuals covariance matrices do not match VAR coefficients matrices');
statsig = m1 == 1;
if statsig
    [Su,cholp] = chol(Su,'lower');
    assert(cholp == 0,'covariance matrix not positive-definite');
else
    assert(m1 == m,'residuals covariance matrices do not match VAR coefficients matrices');
end
Y=zeros(R,n,m);
% Y = zeros(n,m,R);
for r = 1:R % for each realization
    [Y(r,:,:),NN] = genvar_nonstat(A,Su,n,p,m,statsig);
end

function [Y,NN] = genvar_nonstat(A,Su,n,p,m,statsig)

% initialise to Gaussian white noise

if statsig
    Y = mvnrnd(zeros(1,n),Su,m)'; % "Su" is actually Cholesky matrix
else
    Y = zeros(n,m);
    for t = 1:m % for each time step
        [C,cholp] = chol(Su(:,:,t),'lower');
        assert(cholp == 0,'covariance matrix not positive-definite at time step %d',t);
        Y(:,t) = C*randn(n,1);
    end
end
NN=Y;
% loop through time steps
for t = p+1:m   % for each time step
    for k = 1:p % for each lag
        Y(:,t) = Y(:,t) + A(:,k,t)*Y(:,t-k);
    end
end
