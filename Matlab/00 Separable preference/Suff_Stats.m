%%   MPCs
[meshes.m,meshes.h] = ndgrid(grid.m,grid.h);

% NW=par.gamma/(1+par.gamma).*(par.N/par.H).*par.W;
NW=(L_new/par.H).*par.W;
WW=NW.*ones(mpar.nm,mpar.nh); %Wages
WW(:,end)=par.PROFITS*par.profitshare;

% MPC
WW_h=squeeze(WW(1,:)); WW_h_mesh=squeeze(WW(:,:).*meshes.h);
grid_h_aux=grid.h;
MPC_m = zeros(mpar.nm,mpar.nh);

for hh=1:mpar.nh
%     MPC_m(:,hh)=gradient(squeeze(c_guess(:,hh)))./gradient(grid.m)';
    MPC_m(:,hh)=gradient(squeeze(c_new(:,hh)))./gradient(grid.m)'; % MPC_m_ is same with MPC_m
end


MPC_h = zeros(mpar.nm,mpar.nh);

for mm=1:mpar.nm
%     MPC_h(mm,:)=gradient(log(c_guess(mm,:)))./gradient(log(WW_h.*grid_h_aux));
   MPC_h(mm,:)=gradient(c_new(mm,:))./gradient(WW_h.*grid_h_aux);
end

NNP = meshes.m;
URE = inc.labor + meshes.m - c_new;

MPC_m=min(MPC_m,1); % prevent MPC_a larger than 1
%% Sufficient statistics in Auclert

% 1. Aggregate income channel: EI[Yi/Y MPC_m]

Inc_wt_MPC = sum(sum( WW_h_mesh/sum(sum(joint_distr.*WW_h_mesh)).*MPC_m.*joint_distr ));
% income includes interest rate income?

% 3. Fisher channel
Redist_elas_P = sum(sum(MPC_m.*NNP.*joint_distr)) - sum(sum(MPC_m.*joint_distr))*sum(sum(NNP.*joint_distr));

% 4. Interest rate exposure channel
Redist_elas_R = sum(sum(MPC_m.*URE.*joint_distr)) - sum(sum(MPC_m.*joint_distr))*sum(sum(URE.*joint_distr));

% 5. substitution channel
sig_i = par.xi^(-1)*c_new./c_new; % c_guess: composite goods consumption, c_new: final goods consumption

Hick_scaling = sum(sum(sig_i.*(1-MPC_m).*c_new.*joint_distr));

