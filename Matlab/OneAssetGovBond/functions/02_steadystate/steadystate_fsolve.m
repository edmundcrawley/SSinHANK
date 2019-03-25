function [c_guess,m_star,joint_distr_new,W_fc,Profits_fc,Output,N,grid,excess,par,C_agg,C_ind,inc,X_agg]...
    =steadystate_fsolve(P_H,grid,mpar,par,meshes)
   
    
fex=@(Guess)excessB(Guess,grid,P_H,mpar,par,meshes);

Guess=[par.RB];

options=optimset('Display','off','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',4000000);
  
SS_Returns= fsolve(fex,Guess,options);
  
par.RB=SS_Returns;

  %% Update
  [excess,c_guess,m_star,joint_distr_new,W_fc,Profits_fc,Output,N,par,grid,inc]...
                          =excessB(SS_Returns,grid,P_H,mpar,par,meshes);
  
  grid.N=N;
  aux_x = par.tau*N*W_fc*meshes.h/par.H/(1+par.gamma);
  aux_x(:,end)=0;
  C_ind=c_guess+aux_x;
  C_agg = joint_distr_new(:)'*C_ind(:);
  X_agg = joint_distr_new(:)'*c_guess(:);
end

function [excess,c_guess,m_star,joint_distr,W_fc,Profits_fc,Output,N,par,grid,inc]...
    = excessB(Guess, grid,P_H,mpar,par,meshes)

par.RB=Guess(1);

K=grid.K;

    mc   =  par.mu - (par.beta * log(par.PI) - log(par.PI))/par.kappa;
    N    =  (par.tau.*par.alpha.*K^(1-par.alpha).*mc).^(1/(1-par.alpha+par.gamma));   
    W_fc =  par.alpha .*mc.* (K./N).^(1-par.alpha);        
    
    Y = (N).^(par.alpha).*K.^(1-par.alpha);
    Output = Y;
    
    Profits_fc = 1*(1-mc)*Y - Y.*(1/(1-par.mu))./par.kappa./2 .*log(par.PI).^2; 
    
    NW=par.gamma/(1+par.gamma).*(N/par.H).*W_fc;
    WW=NW*ones(mpar.nm,mpar.nh); %Wages
    WW(:,end)=Profits_fc*par.profitshare;
    RBRB = (par.RB+(meshes.m<0)*par.borrwedge)./par.PI;
    
    [c_guess,inc]=policyguess(meshes,WW,RBRB,par,grid,N,W_fc,P_H,Profits_fc,Output);
    count=0;
  
    % 5) Solve Policies and Joint Distribution
    disp('Solving household problem by EGM')
    tic
%% policies_SS
    P=P_H;
    money_expense  = repmat(grid.m',[1 mpar.nh]); % money_expense is new saving
    distC = 99999;
  
    count=0;
    while max([distC])>mpar.crit
        count=count+1;

        % Step 1: Update policies for only money adjustment
        mutil_c = 1./(c_guess.^par.xi); % marginal utility at consumption policy no adjustment

        mutil_c=RBRB.*mutil_c; %take return on money into account
        aux=reshape(permute(mutil_c,[2 1]),[mpar.nh mpar.nm]);
        % form expectations
        EMU_aux = par.beta*permute(reshape(P*aux,[mpar.nh mpar.nm]),[2 1]);

        c_aux = 1./(EMU_aux.^(1/par.xi)); % c_aux is an initial guess for current consumption for a given future consumption c_guess

        % Take borrowing constraint into account
        [c_new,m_star]=EGM_Step1_a(grid,inc,money_expense,c_aux,mpar,par);
        m_star(m_star>grid.m(end)) = grid.m(end);
   
        % Step 6: Check convergence of policies
        distC = max((abs(c_guess(:)-c_new(:))));
    
        % Update c policy guesses
        c_guess=c_new;

    end
    distPOL=[distC];
    toc
    
  disp(([distPOL]));

  disp('Calc Joint Distr')

  [joint_distr]=JDiteration(m_star,P_H,mpar,grid);
  
  joint_distr=reshape(joint_distr,[mpar.nm mpar.nh]);

  aux_x = par.tau*N/par.H*W_fc*meshes.h/(1+par.gamma);
  aux_x(:,end)=0;
  C_ind=c_guess+aux_x;
  C_agg = joint_distr(:)'*C_ind(:);  
  
  AggregateSaving = m_star(:)'*joint_distr(:);
  grid.B=par.BtoY*Output; 
  excess(1)=(grid.B - AggregateSaving);

  
end


