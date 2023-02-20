function [Out_mvDE, npdf]=T_SmvDE_NCDF(X,m,c,stau,designated_signal,threshold)
% This group of scripts calculate stratified multiscale multivariate dispersion entropy
% (SmvMDE) using normal cumulative distribution function (NCDF).

% This is the (1) script used for single scale and the 1 iteration of multiscale operations
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

%% Step 1: Normalization using NCDF and Formulation of "classified" series
CH=size(X,1); % number of channels
N=size(X,2); % length of each signal

for i_CH=1:CH
    x=X(i_CH,:);
    % time series is normalized between (0,1)
    % Normal cumulative distribution function (NCDF) is used 
    % to map x into y from 0 to 1 
    sigma=std(x); 
    mu=mean(x);
    x=normcdf(x,mu,sigma); % Normalization occurs at this line
    
    x(x==1)=1-1e-10; % Soften edges to avoid computation errors
    x(x==0)=1e-10; % Soften edges
    
    X(i_CH,:)=round((x*c)+0.5); % production of classified time-series
    
    % Reason for adding +0.5
    % Since all values of x are in range of 0 < to < 1 the +0.5 is added to
    % ensure that classes are in correct range for example: 1 to 6 instead of 0 to 5

end


%% Step 2: Formulation of multivariate embedded vectors

tau=stau*ones(1,CH); % time delay per channel - in this version the time
% delay is common accross all channels

M=m*ones(1,CH); % embedding dimension value for each channel - common as above

S_M=sum(M); % calculation of m * number of channels

y=[]; % empty array to store the multivariate embedded vectors 

for j_embd=1:CH % we scan through each channel
    for i_embd=1:N-m+1 % total number of multivariate embedded vectors
        
        % (A) Initially, partial embedded vectors are produced per channel
        temp1(i_embd,:)=X(j_embd,i_embd:tau(j_embd):i_embd+M(j_embd)-1);
        % Explanation: formulate part of the embedding vector starting from the current
        % embedding index, use a time delay step associated with each channel
        % build up to the number of embedded samples per channel
    end
    % (B) combine each partial embedded vector from each channel
    % to formulate a complete  multivariate embedding vector
    
    y=horzcat(y,temp1); % combine each partial embedded vector from each channel
    
    
    temp1=[]; % reset the partial embedding vector variable
end

%% Step 3: Generation of all possible DisEn patterns

%  Generate all possible dispersion patterns using a recursive process:

all_pattern=(1:c)'; % (A) starting with embedding dimension 1

for f=2:m % building up recursively to design all possible combinations
    temp=all_pattern;
    all_pattern=[];
    j=1;
    for w=1:c % (B) formulate batches
        [a,~]=size(temp);
        all_pattern(j:j+a-1,:)=[temp,w*ones(a,1)];
        j=j+a;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pattern_num=c^m; % obtain the number of all possible UNIQUE DisEn patterns

key=zeros(1,pattern_num); % array to store each potential DisEn pattern

% Allocate each pattern to a location in the key array
for i_p=1:pattern_num
    for ii_p=1:m
        key(i_p)=key(i_p)*10+all_pattern(i_p,ii_p);
    end
end

pdf=zeros(1,pattern_num); %initialize frequency array for each DisEn pattern

%% Step 4: Accessing Combinations and Trimming based on Designated Signal and Threshold

% (A) The trimming operation is specific to the Threshold variation of the
% Algorithm

aa=combnk(1:S_M,m);  % matrix that stores all possible combinations of
% samples for accessing each multivariate embedded vector

[rows_aa , ~] = size(aa); % retrieve dimensions of access combinations
% rows_aa corresponds to all the possible combinations of m samples that
% can be accessed from the CH * m samples of each multivariate embedded
% vector

% (B) Calculate designated samples based on the input designated_signal which
% indicates the index of the input channel that will be prioritised. 

% The embedded patterns are generated in the same order as the input
% time-series channels. As a result the first m samples of each
% multivariate embedded vector correspond to the first channel, the next m
% samples to the second channel etc..
designated_samples  = []; % initialize designated_samples as empty 

for i= 1:length(designated_signal) % scan through input array
    designated_samples = cat(2,designated_samples,...
        ( ((designated_signal(i)-1) * m) + 1):( designated_signal(i)*m) );
end 
% (C) The following loops will Trim aa to ensure that only combinations
% that meet and surpass the designated_samples threshold are used for
% the calculation of the output DisEn value

% Initialise an index for the formulation of the trimmed combinations
k_trimmed = 1; % index for the aa_trimmed

for i_aa = 1:rows_aa
    
    counter = 0; % initialise counter for each row
    for j_ds = 1:length(designated_samples) % check for each designated sample       
        counter = counter + sum(aa(i_aa,:) == designated_samples(j_ds));   
    end
    
    
   
    % If the row of aa passes the threshold check feed it to aa_trimmed
    if (counter >= threshold )
        
        aa_trimmed(k_trimmed,:) = aa(i_aa,:);        
        k_trimmed = k_trimmed + 1;
    end
    
end 

[rows_aa_trimmed, ~] = size(aa_trimmed); % retrieve dimensions


pv=zeros(rows_aa_trimmed,N-m+1); % array to store dispersion patterns
% corresponding to each accesing combination (rows_aa_trimmed)
% for each embedded vector (N-m+1)

%%  Step 5: Allocation of accessed sample combinations to dispersion patterns
for i_pv=1:rows_aa_trimmed 
    for iii_p=m:-1:1
        pv(i_pv,:)=pv(i_pv,:)+y(:,aa_trimmed(i_pv,iii_p))'*10^(iii_p-1);
    end
end

%% Step 6: Calculation of dispersion pattern frequency and DisEn value

for i=1:pattern_num 
    [a ,~]=find(pv==key(i)); 
    pdf(i)=length(a); 
end


npdf=pdf/((N-m+1)*rows_aa_trimmed); % normalize the frequency array to obtain probability density function.

% Code for alternative calculation of normalized npdf:
% npdf2=pdf/sum(pdf);
p=npdf(npdf~=0);

Out_mvDE=-sum(p .* log(p))/(log(c^m));
% Original Implementation without c^m normalization: Out_mvDE=-sum(p .* log(p));
end
