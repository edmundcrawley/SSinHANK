% Steady state statistics
    clear targets
    targets.ShareBorrower=sum((grid.m<0).*sum(joint_distr,2)');
    targets.B=grid.m*(sum(joint_distr,2));
    targets.B_over_Y=targets.B/Output;
    targets.Y=Output;
    BCaux_M=sum(sum(joint_distr,2),3);
    targets.m_bc=BCaux_M(1,:);
    targets.m_0=BCaux_M(grid.m==0);
    labortax =(1-par.tau).*W_fc.*N +(1-par.tau).*Profits_fc;
    par.gamma1=targets.B*(1-par.RB)+labortax(1);
    par.W=W_fc(1);
    par.PROFITS=Profits_fc(1);
    par.N=N;
    targets.T=labortax;
%     par.G=targets.B*(1-par.RB/par.PI)+targets.T;
    par.G=targets.B*(1-par.RB/par.PI);
    targets.GtoY=par.G/Output;
