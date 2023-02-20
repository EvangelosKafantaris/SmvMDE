function Out_mvMDE=T_SmvMDE(X,m,c,stau,Scale,designated_signal,threshold)
% This group of scripts calculate stratified multiscale multivariate dispersion entropy
% (SmvMDE) using normal cumulative distribution function (NCDF).

% This is the script that calls the scripts (1) and (2) to provide full
% single and multiscale functionality. Use this script to utilize the
% algorithm in your code

% The T-SmvMDE variation includes two new input parameters named:
% designated_signal and threshold

%
% Inputs:
% X: multivariate signal - a matrix of size nvar (the number of channels) x nsamp (the number of sample points for each channel)
% m: scalar embedding value
% c: number of classes (we usullay set c=5, 6, or 7)
% stau: scalar time lag  value (it is usually equal to 1)

% Inputs Specific to this Variation (T-SmvMDE):

% designated_signal: the set of signals (channels) that are allocated to the "core" stratum and prioritised
% threshold: the minimum number of samples from the designated_signal
%            required when accessing each embedding vector 

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

% If the current script is utilized please remember to reference the
% publication that introduced these novel variations:

% [1] E. Kafantaris, T. -Y. M. Lo and J. Escudero, "Stratified Multivariate Multiscale Dispersion Entropy for Physiological Signal Analysis," 
% in IEEE Transactions on Biomedical Engineering, vol. 70, no. 3, pp. 1024-1035, March 2023, doi: 10.1109/TBME.2022.3207582.

% Evangelos Kafantaris and Javier Escudero
% evangelos.kafantaris@ed.ac.uk and javier.escudero@ed.ac.uk

%%
Out_mvMDE=NaN*ones(1,Scale);
% multivariate dispersion entropy at scale 1.
Out_mvMDE(1)=T_SmvDE_multi_NCDF(X,m,c,stau,designated_signal,threshold);


% multivariate dispersion entropy at temporal scales 2 to maximum scale factor.
sigma=std(X,0,2);
mu=mean(X,2);

for j=2:Scale
    Xs = Multi(X,j);
    Out_mvMDE(j)=T_SmvDE_multi_NCDF_ms(Xs,m,c,mu,sigma,stau,designated_signal,threshold);
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