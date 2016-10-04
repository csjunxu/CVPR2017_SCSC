% Main Function of Coupled Dictionary Learning
% Input:
% Alphap,Alphas: Initial sparse coefficient of two domains
% Xp    ,Xs    : Image Data Pairs of two domains
% Dp    ,Ds    : Initial Dictionaries
% Wp    ,Ws    : Initial Projection Matrix
% par          : Parameters 
%
% Output
% Alphap,Alphas: Output sparse coefficient of two domains
% Dp    ,Ds    : Output Coupled Dictionaries
% Up    ,Us    : Output Projection Matrix for Alpha
% 

function [Alphac, Alphan, Uc, Un, Pn, f] = SCSC_DAP_Learning(Alphac, Alphan, Xc, Xn, Dc, Dn, par)

%% parameter setting
param.lambda        = 	    par.lambda1; % not more than 20 non-zeros coefficients
param.lambda2       =       par.lambda2;
param.mode          = 	    2;       % penalized formulation
param.approx=0;
param.K = par.K;
param.L = par.L;
f = 0;

%% Initialize Us, Up as I

Un = eye(size(Dn, 2));
Uc = eye(size(Dc, 2));

%% Initialize Ps as 0 matrix
Pn = zeros(size(Xn));

%% Iteratively solve A U P

for t = 1 : par.nIter
    
    %% Updating Alphas and Alphap
    f_prev = f;
    Alphan = mexLasso([Xn - Pn; par.sqrtmu * Uc * full(Alphac)], [Dn; par.sqrtmu * Un],param);
    Alphac = mexLasso([Xc; par.sqrtmu * Un * full(Alphan)], [Dc; par.sqrtmu * Uc],param);
   
    %% Updating Ps
    Pn = (Xn - Dn * Alphan) / (1 + par.nup);

    %% Updating Us and Up
    Un = (1 - par.rho) * Un  + par.rho * Uc * Alphac * Alphan' / ( Alphan * Alphan' + par.nu * eye(size(Alphan, 1)));
    Uc = (1 - par.rho) * Uc  + par.rho * Un * Alphan * Alphac' / ( Alphac * Alphac' + par.nu * eye(size(Alphac, 1)));

    %% Find if converge (NEED MODIFICATION)
    P1 = Xc - Dc * Alphac;
    P1 = P1(:)'*P1(:) / 2;
    P2 = par.lambda1 *  norm(Alphac, 1);    
    P3 = Un * Alphan - Uc * Alphac; 
    P3 = P3(:)'*P3(:) / 2;
    P4 = par.nu * norm(Uc, 'fro');
    fp = 1 / 2 * P1 + P2 + par.mu * (P3 + P4);
    
    P1 = Xn - Dn * Alphan - Pn;
    P1 = P1(:)'*P1(:) / 2;
    P2 = par.lambda1 *  norm(Alphan, 1);    
    P3 = Un * Alphan - Uc * Alphac; 
    P3 = P3(:)'*P3(:) / 2;
    P4 = par.nu * norm(Un, 'fro');  
    P5 = par.nup * norm(Pn, 'fro'); 
    fs = 1 / 2 * P1 + P2 + par.mu * (P3 + P4 + P5);
    
    f = fp + fs;
	
    %% if converge then break
    if (abs(f_prev - f) / f < par.epsilon)
        break;
    end
    fprintf('Energy: %d\n',f);
end
    
