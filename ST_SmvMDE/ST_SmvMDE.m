function Out_mvMDE=ST_SmvMDE(X,m,c,stau,Scale,designated_signal,threshold,reduced_weight)
% This group of scripts calculate stratified multiscale multivariate dispersion entropy
% (SmvMDE) using normal cumulative distribution function (NCDF).

% This is the script that calls the scripts (1) and (2) to provide full
% single and multiscale functionality. Use this script to call utilize the
% algorithm in your code

% The ST-SmvMDE variation includes three new input parameters named:
% designated_signal, threshold and reduced_weight

%
% Inputs:
% X: multivariate signal - a matrix of size nvar (the number of channels) x nsamp (the number of sample points for each channel)
% m: scalar embedding value
% c: number of classes (we usullay set c=5, 6, or 7)
% stau: scalar time lag  value (it is usually equal to 1)

% Inputs Specific to this Variation (ST-SmvMDE):

% designated_signal: the set of signals (channels) that are allocated to the "core" stratum and prioritised

% threshold: the minimum number of samples from the designated_signal
%            required when accessing each embedding vector 

% reduced_weight: the weight factor used limit the contribution of
% accessing combinations that do not meet the set threshold. Important, its
% set value should be between 0 and 1 for the algorithm to operate as
% designed.

% Output:
% Out_mvMDE: a scaler value - the SmvDE of X
%

% This variation is based on the original Multivariate Multiscale
% Dispersion Entropy Algorithm (mvMDE) for which the respective references
% are:

% [1] H. Azami, A. Fernandez, and J. Escudero, "Multivariate Multiscale Dispersion Entropy of Biomedical Times Series", Entropy, 2019.
% [2] H. Azami, M. Rostaghi, D. Abasolo, and J. Escudero,"Refined composite multiscale dispersion entropy and its
% application to biomedical signals", IEEE Transactions on Biomedical Engineering, 64(12), 2872-2879, 2017.
%
% Hamed Azami and Javier Escudero
% hmd.azami@gmail.com and javier.escudero@ed.ac.uk
%
%  10-Sep-19

% If the current script is utilized please remember to reference the
% publication that introduced these novel variations:

% [1]


% Evangelos Kafantaris and Javier Escudero
% evangelos.kafantaris@ed.ac.uk and javier.escudero@ed.ac.uk
%  6-April-2021
%%
Out_mvMDE=NaN*ones(1,Scale);
% multivariate dispersion entropy at scale 1.
Out_mvMDE(1)=ST_SmvDE_multi_NCDF(X,m,c,stau,designated_signal,threshold,reduced_weight);


% multivariate dispersion entropy at temporal scales 2 to maximum scale factor.
sigma=std(X,0,2);
mu=mean(X,2);

for j=2:Scale
    Xs = Multi(X,j);
    Out_mvMDE(j)=ST_SmvDE_multi_NCDF_ms(Xs,m,c,mu,sigma,stau,designated_signal,threshold,reduced_weight);
end

%% Define Function for the Coarse-Graining of Input Samples
function M_Data = Multi(Data,S)

%  generate the consecutive coarse-grained time series
%  Input:   Data: time series;
%           S: scale factor

% Output:
%           M_Data: coarse-grained time series at scale factor S

L = size(Data,2);
J = fix(L/S);

for j=1:size(Data,1)
    for i=1:J
        M_Data(j,i) = mean(Data(j,(i-1)*S+1:i*S));
    end
end