def SolveSteadyState(self):

       
        ## Set grid h
        grid = grid
        resultStVar=StochasticsVariance(par, mpar, grid)
       
        P_H = resultStVar['P_H'].copy()
        grid = resultStVar['grid'].copy()
        par = resultStVar['par'].copy()


        regrid = MakeGrid2(mpar, grid, grid['m_min'], grid['m_max'])
           
        grid['m'] = regrid['m'].copy()
        
        meshesm, meshesh =  np.meshgrid(grid['m'],grid['h'],indexing='ij')
        meshes = {'m': meshesm, 'h': meshesh}
        
        fex = lambda Guess : excessB(Guess,grid,P_H,mpar,par,meshes)
        
        Guess = par['RB']
        SS_Return = fsolve(fex,Guess)
        par['RB']=SS_Return

        Results_excessB = excessB_update(SS_Return,grid,P_H,mpar,par,meshes)

        excess = Results_excessB['excess']
        c_new = Results_excessB['c_policy']
        m_star =  Results_excessB['m_policy']
        joint_distr =  Results_excessB['joint_distr']
        W_fc =  Results_excessB['W_fc']
        Profits_fc =  Results_excessB['Profits_fc']
        Output =  Results_excessB['Output']
        L_new =  Results_excessB['L_new']
        grid = Results_excessB['grid']
        inc = Results_excessB['inc']

        grid['N']=L_new.flatten('F').transpose().dot(joint_distr.flatten('F').transpose())
        
        C_agg = np.sum(np.sum(np.multiply(joint_distr,c_new)))

        
        ## SS_stats
        tgSB = np.sum((grid['m']<0)*np.sum(joint_distr.copy(),1))
        tgB = grid['m'].copy().dot(np.sum(joint_distr.copy(),1))
        tgBY = tgB/Output
        BCaux_M = np.sum(joint_distr,1) # 
        tgm_bc = BCaux_M[0,:]
        tgm_0 = BCaux_M[grid['m']==0]
        
        tax = (1.-par['tau'])*W_fc*N +(1.-par['tau'])*Profits_fc
        tgG=tgB*(1-par['RB'])+tax
        tgW=W_fc
        tgPROFITS=Profits_fc
        tgN = N
        tgGtoY = tgG/Output
        tgT = tax
        
              
        # Net wealth Gini
        Ginipdf = np.array(np.sum(joint_distr.copy(),1).transpose())
        Ginicdf = np.array(np.cumsum(joint_distr.copy(),axis=0))
        
        grid_m_aux = abs(min(grid['m'].copy()))+grid['m'].copy()+0.1                
        #GiniS_aux0 = np.cumsum(Ginipdf.copy()*grid['m'].copy())
        GiniS_aux = np.cumsum(Ginipdf.copy()*grid_m_aux.copy())       
        
        GiniS = np.concatenate((np.array([0.]),GiniS_aux.copy()),axis=0)
        tgGini = 1.-np.sum(Ginipdf*(GiniS.copy()[:-1]+GiniS.copy()[1:]))/GiniS.copy()[-1]
        
        targets = {'ShareBorrower' : tgSB,
                   'B': tgB,
                   'BY': tgBY,
                   'Y': Output,
                   'm_bc': tgm_bc,
                   'm_0': tgm_0,
                   'G': tgG,
                   'W': tgW,
                   'PROFITS': tgPROFITS,
                   'N': tgN,
                   'GtoY': tgGtoY,
                   'T': tgT,
                   'Gini': tgGini
                   }   
        
        
        ## Prepare state and controls        
        grid['B'] = np.sum(grid['m']*np.sum(joint_distr,1))
        
        # Calculate Marginal Values of Assets (m)
        RBRB = (par['RB']+(meshes['m'].copy()<0)*par['borrwedge'])/par['PI']
        
        # Liquid Asset
        mutil_c = 1./(c_guess.copy()**par['xi'])
        Vm = RBRB.copy()*mutil_c.copy()
        Vm = np.reshape(Vm.copy(),(mpar['nm'],mpar['nh']))
        
        
        ## Produce non-parametric Copula
        cum_dist = np.cumsum(np.cumsum(joint_distr,axis=0),axis=1)
        marginal_m = np.cumsum(np.squeeze(np.sum(joint_distr,axis=1)))
        marginal_h = np.cumsum(np.squeeze(np.sum(joint_distr,axis=0)))
        
        
        Copula = interp2d(marginal_m.copy(), marginal_h.copy(), cum_dist.copy().transpose(),kind='cubic')
        
       
        return {'par':par,
                'mpar':mpar,
                'grid':grid,
                'Output':Output,
                'targets':targets,
                'Vm': Vm,
                'joint_distr': joint_distr,
                'Copula': Copula,
                'c_policy': c_guess,
                'm_policy': m_star,
                'mutil_c': mutil_c,
                'P_H' : P_H,
                'inc' : inc,
                'Profits_fc' : Profits_fc,
                'C_agg' : C_agg,
                'C_ind': C_ind,
                'X_agg' : X_agg
                }


def FactorReturns(meshes, grid, par, mpar,L_guess):

        ## GHH preferences
        mc = par['mu'] - (par['beta'] * np.log(par['PI']) - np.log(par['PI']))/par['kappa']

#        L_guess = np.ones((mpar['nm'],mpar['nh']))
#        L_guess[:,-1] = 0
        
#        w =  par['alpha'] *mc * (grid['K']/N) **(1-par['alpha'])
        w = mc
        Y = np.sum(np.sum(L_guess))/np.sum(np.sum(np.ones((mpar['nm'],mpar['nh']))))
        Profits_fc = (1-mc)*Y
        
#        Y = N**par['alpha']*grid['K']**(1-par['alpha'])
    
        WW = w*L_guess
        WW[:,-1] = Profits_fc * par['profitshare']
        RBRB = (par['RB']+(meshes['m']<0)*par['borrwedge'])/par['PI']
    
        return {'L_guess':L_guess, 'W_fc':w, 'Profits_fc':Profits_fc,'WW':WW,'RBRB':RBRB,'Y':Y}
    

def PolicyGuess(meshes, WW, RBRB, par, mpar,grid, W_fc, P_H, Profits_fc,Output):

        
        jd_aux= np.linalg.matrix_power(P_H.copy(),1000)
        jd_aux=jd_aux[0,:]

        inclabor = par['tau']*WW*meshes['h']/par['H'].copy()
        incmoney = RBRB*meshes['m'].copy()
        #incprofits = sum((1-par['tau'])*par['gamma']/(1+par['gamma'])*(N/par['H'])*W_fc*grid['h'][0:-1]*jd_aux[0:-1]) + (1-par['tau'])*Profits_fc*par['profitshare']*jd_aux[-1]
          
#        incprofits = sum((1-par['tau'])*(N/par['H'])*W_fc*grid['h'][0:-1]*jd_aux[0:-1]) \
#                         + (1-par['tau'])*Profits_fc*par['profitshare']*jd_aux[-1] \
#                         - (RBRB-1.0)*par['BtoY']*Output

           
#        inc = {'labor': inclabor.copy(),
#               'money': incmoney.copy(),
#               'profits': incprofits}
        
        inc = {'labor': inclabor.copy(),'money': incmoney.copy()}        
    
#        c_guess = inc['labor'].copy()+np.maximum(inc['money'],0.) + inc['profits']
        c_guess = inc['labor'].copy()+np.maximum(inc['money'],0.)        
    
        return {'c_guess':c_guess, 'inc': inc}       


def PoliciesSS(c_guess, grid, inc, RBRB, P, mpar, par, meshes, Profits_fc, W_fc):

             
        ## Apply EGM to slove for optimal policies and marginal utilities
        money_expense = np.transpose(np.tile(grid['m'],(mpar['nh'],1)))
        distC = 99999.
    
        count = 0
    
        while np.max((distC)) > mpar['crit']:
        
               count = count+1
        
               ## update policies
               mutil_c = 1./(c_guess.copy()**par['xi']) # marginal utility at consumption policy
        
               aux = np.reshape(np.ndarray.transpose(mutil_c.copy(),(1,0)),(mpar['nh'],mpar['nm']) )
        
               # form expectations
               EMU_aux = par['beta']*RBRB*np.ndarray.transpose(np.reshape(P.copy().dot(aux.copy()),(mpar['nh'],mpar['nm'])),(1,0))
        
               c_aux = 1./(EMU_aux.copy()**(1/par['xi']))
        
               # Take budget constraint into account
               resultEGM = EGM(grid, inc, money_expense, c_aux, mpar, par, meshes, Profits_fc, W_fc)
               c_new = resultEGM['c_update'].copy()
               m_star = resultEGM['m_update'].copy()
        
               m_star[m_star>grid['m'][-1]] = grid['m'][-1] # no extrapolation
        
               ## check convergence of policies
               distC = np.max((np.abs(c_guess.copy()-c_new.copy())))
        
               # update c policy guesses
               c_guess = c_new.copy()
        
        distPOL = distC

        return {'c_new':c_new, 'm_star':m_star, 'distPOL':distPOL}     


def EGM(grid, inc, money_expense, c_aux, mpar, par, meshes, Profits_fc, W_fc):

    
        ## EGM: Calculate assets consistent with choices being (m')
        # Calculate initial money position from the budget constraint,
        # that leads to the optimal consumption choice
        mmin = grid['m'][0]
#        mmax = grid['m'][-1]
        
        cpbind = ((par['RB'])*meshes['m']+mmin)/2+np.sqrt(((par['RB'])*meshes['m']+mmin)**(2)+4*(meshes['h']/par['H']*W_fc)**(1+1/par['gamma']))/2

        cpbind[:,-1] = ((par['RB'])*meshes['m'][:,-1]+mmin+Profits_fc*par['profitshare']);
                  
                 
#        m_star = c_aux + money_expense - inc['labor'] -inc['profits']
        m_star = c_aux + money_expense - inc['labor']       
        
        RR = (par['RB']+(m_star.copy()<0.)*par['borrwedge'])/par['PI']
        m_star = m_star.copy()/RR
    
        # Identify binding constraints
        binding_constraints = (money_expense < np.tile(m_star[0,:],(mpar['nm'], 1)))

        # Consumption when drawing assets m' to zero: Eat all Resources
#        Resource = inc['labor']  + inc['money'] + inc['profits']
    
        ## Next step : interpolate on grid
        c_update = np.zeros((mpar['nm'],mpar['nh']))
        m_update = np.zeros((mpar['nm'],mpar['nh']))
        
        for hh in range(mpar['nh']):
            
            Savings = interp1d(m_star[:,hh].copy(), grid['m'], fill_value='extrapolate')
            m_update[:,hh] = Savings(grid['m'])
            Consumption = interp1d(m_star[:,hh].copy(), c_aux[:,hh], fill_value='extrapolate')
            c_update[:,hh] = Consumption(grid['m'])
        
        
#        binding_constraints = c_update > Resource
        c_update[binding_constraints] = cpbind[binding_constraints]
        m_update[binding_constraints] = np.min(grid['m'])
        
        #c_update[binding_constraints] = Resource[binding_constraints]-grid['m'][0]
        #m_update[binding_constraints] = np.min((grid['m']))            

        return {'c_update': c_update, 'm_update': m_update}
       
        
def excessB(Guess, grid,P_H,mpar,par,meshes):
        
        par['RB']=Guess
        
        grid['K']=1
        
        mc = par['mu'] - (par['beta'] * np.log(par['PI']) - np.log(par['PI']))/par['kappa']
        
        L_guess = np.ones((mpar['nm'],mpar['nh']))
        L_guess[:,-1] = 0
        
        ################################################################################################
        count_L = 0
        dist_L = 99999.
        
        while np.max((dist_L)) > mpar['crit']:

             resultFactReturn = FactorReturns(meshes, grid, par, mpar,L_guess)
             L_guess = resultFactReturn['L_guess']
             W_fc = resultFactReturn['W_fc']
             Profits_fc =resultFactReturn['Profits_fc']
             WW = resultFactReturn['WW'].copy()
             RBRB = resultFactReturn['RBRB'].copy()
             Output = resultFactReturn['Y']
            
             resultPolGuess = PolicyGuess(meshes,WW,RBRB,par,mpar,grid,W_fc,P_H,Profits_fc,Output)
             c_guess = resultPolGuess['c_guess'].copy()
             inc = resultPolGuess['inc'].copy()


             print('Solving household problem by EGM')
             
             start_time = time.clock()
           
             resultPolicySS = PoliciesSS(c_guess, grid, inc, RBRB, P_H, mpar, par, meshes, Profits_fc, W_fc)
             c_new = resultPolicySS['c_new'].copy()
             m_star = resultPolicySS['m_star'].copy()
             distC = resultPolicySS['distPOL']
        
             end_time = time.clock()
             print('Elapsed time is ',  (end_time-start_time), ' seconds.')
             print(distC)
             
             meshes_E = np.hstack((meshes['h'][:,0:-1],np.zeros((mpar['nm'],1))))
             L_new = (np.multiply(c_new**(-par['xi']),meshes_E/par['H'])*W_fc)**(1/par['gamma'])
             dist_L = np.max((np.abs(L_guess.copy()-L_new.copy())))
             
             print(dist_L)
 
             L_guess = L_new
             
             count_L=count_L+1

             print('Calc Joint Distr')
             joint_distr = JDiteration(m_star, P_H, mpar, grid)
             
             N = L_new.flatten('F').transpose().dot(joint_distr.flatten('F').transpose())
             Output = N
             Profits_fc = (1-mc)*Output
        #################################################################################################        
        
        print('Two loops for C and L, done.')

        Output = np.squeeze(np.asarray(Output))
        joint_distr = np.reshape(joint_distr.copy(),(mpar['nm'],mpar['nh']),order='F')
        AggregateSavings = m_star.flatten('F').transpose().dot(joint_distr.flatten('F').transpose())
        AggregateSavings = np.asarray(AggregateSavings).reshape(-1)
        ExcessA = AggregateSavings-par['BtoY']*Output      
        
        return ExcessA
        

def excessB_update(SS_return, grid,P_H,mpar,par,meshes):
        
        par['RB']=SS_return
        
        grid['K']=1
        
        mc = par['mu'] - (par['beta'] * np.log(par['PI']) - np.log(par['PI']))/par['kappa']        
        
        L_guess = np.ones((mpar['nm'],mpar['nh']))
        L_guess[:,-1] = 0
        
        ################################################################################################
        count_L = 0
        dist_L = 99999.
        
        while np.max((dist_L)) > mpar['crit']:

             resultFactReturn = FactorReturns(meshes, grid, par, mpar,L_guess)
             L_guess = resultFactReturn['L_guess']
             W_fc = resultFactReturn['W_fc']
             Profits_fc =resultFactReturn['Profits_fc']
             WW = resultFactReturn['WW'].copy()
             RBRB = resultFactReturn['RBRB'].copy()
             Output = resultFactReturn['Y']            
            
             resultPolGuess = PolicyGuess(meshes,WW,RBRB,par,mpar,grid,W_fc,P_H,Profits_fc,Output)
             c_guess = resultPolGuess['c_guess'].copy()
             inc = resultPolGuess['inc'].copy()


             print('Solving household problem by EGM')
             
             start_time = time.clock()
           
             resultPolicySS = PoliciesSS(c_guess, grid, inc, RBRB, P_H, mpar, par, meshes, Profits_fc, W_fc)
             c_new = resultPolicySS['c_new'].copy()
             m_star = resultPolicySS['m_star'].copy()
             distC = resultPolicySS['distPOL']
        
             end_time = time.clock()
             print('Elapsed time is ',  (end_time-start_time), ' seconds.')
             print(distC)
             
             meshes_E = np.hstack((meshes['h'][:,0:-1],np.zeros((mpar['nm'],1))))
             L_new = (np.multiply(c_new**(-par['xi']),meshes_E/par['H'])*W_fc)**(1/par['gamma'])
             dist_L = np.max((np.abs(L_guess.copy()-L_new.copy())))
             
             print(dist_L)
 
             L_guess = L_new
             
             count_L=count_L+1

             print('Calc Joint Distr')
             joint_distr = JDiteration(m_star, P_H, mpar, grid)
             
             N = L_new.flatten('F').transpose().dot(joint_distr.flatten('F').transpose())
             Output = N
             Profits_fc = (1-mc)*Output
        #################################################################################################        
        
        print('Two loops for C and L, done.')   
        
        joint_distr = np.reshape(joint_distr.copy(),(mpar['nm'],mpar['nh']),order='F')
        AggregateSavings = m_star.flatten('F').transpose().dot(joint_distr.flatten('F').transpose())
        AggregateSavings = np.asarray(AggregateSavings).reshape(-1)
        ExcessA = AggregateSavings-par['BtoY']*Output          
        
        return{'excess':ExcessA,'c_policy':c_new,'m_policy':m_star,'joint_distr':joint_distr,
               'W_fc':W_fc,'Profits_fc':Profits_fc,'Output':Output,'L_new':L_new,'par':par,'grid':grid,'inc':inc}

    
def JDiteration(m_star, P_H, mpar, grid):
        '''
        Iterates the joint distribution over m,k,h using a transition matrix
        obtained from the house distributing the households optimal choices. 
        It distributes off-grid policies to the nearest on grid values.
        
        parameters
        ------------
        m_star :np.array
            optimal m func
        P_H : np.array
            transition probability    
        mpar : dict
             parameters    
        grid : dict
             grids
             
        returns
        ------------
        joint_distr : np.array
            joint distribution of m and h
        
        '''
        ## find next smallest on-grid value for money and capital choices
        weight11  = np.zeros((mpar['nm'], mpar['nh'],mpar['nh']))
        weight12  = np.zeros((mpar['nm'], mpar['nh'],mpar['nh']))
    
        # Adjustment case
        resultGW = GenWeight(m_star, grid['m'])
        Dist_m = resultGW['weight'].copy()
        idm = resultGW['index'].copy()
        
        idm = np.transpose(np.tile(idm.flatten(order='F'),(mpar['nh'],1)))
        idh = np.kron(range(mpar['nh']),np.ones((1,mpar['nm']*mpar['nh'])))
        idm = idm.copy().astype(int)
        idh = idh.copy().astype(int)
        
        
        index11 = np.ravel_multi_index([idm.flatten(order='F'), idh.flatten(order='F')],(mpar['nm'],mpar['nh']),order='F')
        index12 = np.ravel_multi_index([idm.flatten(order='F')+1, idh.flatten(order='F')],(mpar['nm'],mpar['nh']),order='F')
        
        
        for hh in range(mpar['nh']):
        
            # Corresponding weights
            weight11_aux = (1.-Dist_m[:,hh].copy())
            weight12_aux =  (Dist_m[:,hh].copy())
    
            # Dimensions (mxk,h',h)   
            weight11[:,:,hh]=np.outer(weight11_aux.flatten(order='F'),P_H[hh,:].copy())
            weight12[:,:,hh]=np.outer(weight12_aux.flatten(order='F'),P_H[hh,:].copy())
        
            
        
        weight11 = np.ndarray.transpose(weight11.copy(),(0,2,1))
        weight12 = np.ndarray.transpose(weight12.copy(),(0,2,1))
        
        rowindex = np.tile(range(mpar['nm']*mpar['nh']),(1,2*mpar['nh']))
        
        
        H = sp.coo_matrix((np.concatenate((weight11.flatten(order='F'),weight12.flatten(order='F'))), 
                       (rowindex.flatten(), np.concatenate((index11.flatten(order='F'),index12.flatten(order='F'))))),shape=(mpar['nm']*mpar['nh'], mpar['nm']*mpar['nh']))
        
        ## Joint transition matrix and transitions
        
        distJD = 9999.
        countJD = 1
    
        eigen, joint_distr = sp.linalg.eigs(H.transpose(), k=1, which='LM')
        joint_distr = joint_distr.copy().real
        joint_distr = joint_distr.copy().transpose()/(joint_distr.copy().sum())
        
        while (distJD > 10**(-14) or countJD<50) and countJD<10000:
        
            joint_distr_next = joint_distr.copy().dot(H.copy().todense())
                                       
            joint_distr_next = joint_distr_next.copy()/joint_distr_next.copy().sum(axis=1)
            distJD = np.max((np.abs(joint_distr_next.copy().flatten()-joint_distr.copy().flatten())))
            
            countJD = countJD +1
            joint_distr = joint_distr_next.copy()
                 
            
        return joint_distr
          

         

def StochasticsVariance(par, mpar, grid):

    
        # First for human capital
        TauchenResult = Tauchen(par['rhoH'], mpar['nh']-1, 1., 0., mpar['tauchen'])
        hgrid = TauchenResult['grid'].copy()
        P_H = TauchenResult['P'].copy()
        boundsH = TauchenResult['bounds'].copy()
    
        # correct long run variance for human capital
        hgrid = hgrid.copy()*par['sigmaH']/np.sqrt(1-par['rhoH']**2)
        hgrid = np.exp(hgrid.copy()) # levels instead of logs
    
        grid['h'] = np.concatenate((hgrid,[1]), axis=0)
        
        P_H = Transition(mpar['nh']-1, par['rhoH'], np.sqrt(1-par['rhoH']**2), boundsH)
    
        # Transitions to enterpreneur state
        P_H = np.concatenate((P_H.copy(),np.tile(mpar['in'],(mpar['nh']-1,1))), axis=1)
        lastrow = np.concatenate((np.tile(0.,(1,mpar['nh']-1)),[[1-mpar['out']]]), axis=1)
        lastrow[0,int(np.ceil(mpar['nh']/2))-1] = mpar['out'] 
        P_H = np.concatenate((P_H.copy(),lastrow), axis=0)
        P_H = P_H.copy()/np.transpose(np.tile(np.sum(P_H.copy(),1),(mpar['nh'],1)))
    
        Paux = np.linalg.matrix_power(P_H.copy(),1000)
        hh = Paux[0,:mpar['nh']-1].copy().dot(grid['h'][:mpar['nh']-1].copy())
    
        par['H'] = hh # Total employment
        par['profitshare'] = Paux[-1,-1]**(-1) # Profit per household
        grid['boundsH'] = boundsH
     
    
        return {'P_H': P_H, 'grid':grid, 'par':par}

