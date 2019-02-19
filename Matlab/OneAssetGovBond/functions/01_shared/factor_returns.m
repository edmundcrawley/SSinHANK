function [ N,w,Profits_fc,WW,RBRB,Y ] = factor_returns(meshes,grid,par,mpar)
%factor_returns
  
%%  GHH Preferences    
  
    mc =  par.mu; % in SS
    N    =  (1*par.tau.*par.alpha.*grid.K.^(1-par.alpha).*mc).^(1/(1-par.alpha+par.gamma));   
    w = 1*par.alpha .*mc.* (grid.K./N).^(1-par.alpha);
     
    Y = 1*(N).^(par.alpha).*grid.K.^(1-par.alpha);
    Profits_fc = (1-mc)*Y;
 %%
    NW=par.gamma/(1+par.gamma).*(N/par.H).*w;
    WW=NW*ones(mpar.nm,mpar.nh); %Wages
    WW(:,end)=Profits_fc*par.profitshare;
    RBRB = (par.RB+(meshes.m<0)*par.borrwedge)./par.PI;


end

