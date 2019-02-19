function [P_H,grid,par]=ExTransitions(S,grid,mpar,par)

% Generate transition probabilities and grid
% [hgrid,P_H,boundsH] = Tauchen(par.rhoH,mpar.nh-1); % LR variance = 1
% [hgrid,P_H,boundsH] = Tauchen(par.rhoH,mpar.nh,1,0,'equi'); % LR variance = 1
% Correct long run variance for *human capital*
% hgrid               = hgrid*par.sigmaH/sqrt(1-par.rhoH^2);
% hgrid               = exp(hgrid); % Levels instead of Logs


aux   = sqrt(S)*sqrt(1-par.rhoH^2); % Short run std

P     = transition(mpar.nh-1,par.rhoH,aux,grid.boundsH);

% grid.h=[hgrid 1];


P_H=[P repmat(mpar.in,[mpar.nh-1 1])];
lastrow=[repmat(0,[1, mpar.nh-1]) 1-mpar.out];
lastrow(ceil(mpar.nh/2))=mpar.out;
P_H=[P_H; lastrow];
P_H=P_H./repmat(sum(P_H,2),[1 mpar.nh]);

% Paux=P_H^1000^1000;
% hh=Paux(1,1:mpar.nh-1)*grid.h(1:mpar.nh-1)';
% par.H=hh(1); %Total Employment

end