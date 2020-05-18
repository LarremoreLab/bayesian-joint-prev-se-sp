import numpy as np
import pickle
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri

def evaluate_exp_min_over_grid(grid,pr,se,sp,N,n_samps):
    f = np.inf*np.ones(np.shape(grid[0]))
    steps = np.shape(grid[0])[0]
    for i in range(steps):
        for j in range(steps):
            Npos = grid[0][i,j].astype(int)
            Nneg = grid[1][i,j].astype(int)
            Nfield = N-Npos-Nneg
            if (Nfield <=0) or (Npos<0) or (Nneg<0):
                continue
            f[i,j] = expected_mean_width(Nfield,Npos,Nneg,pr,se,sp,n_samps)
            print(f[i,j],Nfield,Npos,Nneg)
    return grid[0],grid[1],f

def expected_mean_width(Nfield,Npos,Nneg,prev,sens,spec,n_samps):
    p = p_seropositive_r(prev,sens,spec)
    epos = int(np.round(Nfield*p))
    eneg = int(Nfield-epos)
    etn = int(np.round(Nneg*spec))
    efp = int(Nneg-etn)
    etp = int(np.round(Npos*sens))
    efn = int(Npos-etp)
   
    samples = np.array(mcmc(samps=n_samps,
        pos = epos,
        n = Nfield,
        tp = etp,
        tn = etn,
        fp = efp, 
        fn = efn))
    
    previ = samples[:,0]
    CI = np.percentile(previ,[5,95])
    return CI[1]-CI[0]

def p_seropositive_r(r,sensitivity,specificity):
    '''
    p_seropositive_r(r,sensitivity,specificity)
        Computes the probability that a test comes back seropositive, given
        the true seroprevalence in the population and test specs

    r - true seroprevalence
    sensitivity - sensitivity
    specificity - specificity
    '''
    return r*sensitivity+(1-r)*(1-specificity)


import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
ro.r("source('singleSERO_uncertainTEST.R')")
mcmc = ro.globalenv['sample_posterior_r_mcmc_testun']


N=1000
se = 0.93
sp = 0.98
dx = 20
x = np.arange(dx,601,dx)
grid = np.meshgrid(x,x)
n_samps = 10000

pr = 0.15
X,Y,W = evaluate_exp_min_over_grid(grid,pr,se,sp,N,n_samps)
pickle.dump([X,Y,W],open('W2mcmc2.pkl','wb'))
