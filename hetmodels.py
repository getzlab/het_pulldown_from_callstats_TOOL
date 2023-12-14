import numpy as np
import pandas as pd
from scipy.stats import beta, binom, betabinom
from scipy.optimize import minimize

def run_snp_mixture_model(nref, nalt, include_noise_component=True,fix_noise_component=None,alpha=None, nit=50):
    """ Fit SNP mixture model

    nref -- number of reference bases
    nalt -- number of alternate bases
    alpha -- dirichlet prior on het/homref/homalt/other probabilities
    nit -- number of iterations to run EM algorithm
    """
    print('running mixture model')

    if alpha is None:
        f_alt = sum(nalt)/(sum(nref)+sum(nalt))
        alpha = 10000 * np.array([2*f_alt*(1-f_alt),(1-f_alt)**2,f_alt**2])

        if include_noise_component and (fix_noise_component is None):
            alpha = np.hstack([alpha,np.array([.1])])

    # Initialize parameters to be optimized
    pi = alpha / alpha.sum()

    if fix_noise_component is not None:
        pi = np.hstack([pi * (1 - fix_noise_component),
                        np.array([fix_noise_component])])

    skew = .5
    epsilon = 1e-3

    # Perform EM
    for i in range(0, nit):
        # E-step
        snp_post = method_mixture_postprob(nref, nalt, epsilon=epsilon, skew=skew, pi=pi,include_noise_component=include_noise_component)

        # M-step

        if fix_noise_component is not None:
            # If fixing noise frequency, then just update het/hom frequencies
            class_sum = snp_post.sum(axis=1)[0:3] + alpha
            pi = np.hstack([class_sum/class_sum.sum()*(1-fix_noise_component),
                            np.array([fix_noise_component])])
        else:
            # Otherwise, update all of pi (may or may not include noise term
            pi = snp_post.sum(axis=1) + alpha
            pi = pi / pi.sum()

        res = minimize(lambda x: method_mixture_Eloss(snp_post, nref, nalt, x[0], x[1],include_noise_component),
                       x0=(np.log(epsilon), skew),
                       bounds=[(np.log(1e-20), np.log(1e-2)), (.4, .6)])

        epsilon = np.exp(res.x[0])
        skew = res.x[1]

        print(f'skew : {skew}, epsilon : {epsilon}, pi : {pi}')

    colnames = ['prob_het','prob_homref','prob_homalt']
    if include_noise_component:
        colnames = colnames + ['prob_other']

    p_df = pd.DataFrame(snp_post.T,columns = colnames)

    return ({'snp_prob':p_df,'epsilon': epsilon, 'skew': skew, 'pi': pi})

# Returns a 4 x N matrix with the posterior probability each SNP is a het, hom_ref, hom_alt, or none of the above
def method_mixture_postprob(nref,nalt,epsilon,skew,pi,include_noise_component):

    nref=nref.reshape(-1)
    nalt=nalt.reshape(-1)

    lik = method_mixture_likelihood(nref,nalt,epsilon=epsilon,skew=skew,include_noise_component=include_noise_component)
    
    pi_shape = np.ones(len(lik.shape)).astype(int)
    pi_shape[0]=4 if include_noise_component else 3
    pi = pi.reshape(pi_shape)
    
    post_prob = (lik*pi) / (lik * pi).sum(axis=0)

    return(post_prob)

# Returns a 4 x N matrix with the model likelihood each SNP is a het, hom_ref, hom_alt, or none of the above
def method_mixture_likelihood(nref,nalt,epsilon,skew,include_noise_component):
    
    N = nref+nalt
    
    het_lik = binom.pmf(nalt,N,(1-epsilon)*skew + epsilon*(1-skew))
    homref_lik = binom.pmf(nalt,N,epsilon)
    homalt_lik = binom.pmf(nalt,N,1-epsilon)

    lik = np.stack([het_lik,homref_lik,homalt_lik])

    if include_noise_component:
        unif_lik = betabinom.pmf(nalt, N, 1, 3)
        lik = np.vstack([lik,unif_lik.reshape(1,len(nalt))])

    return(lik)
    
# Negative expected log posterior to optimize for the M-step
def method_mixture_Eloss(snp_post,nref,nalt,log_epsilon,skew,include_noise_component,min_lik=1e-50):

    lik = method_mixture_likelihood(nref,nalt,np.exp(log_epsilon),skew,include_noise_component)
    lik = np.maximum(min_lik,lik)
    lik = lik.reshape(snp_post.shape[0], lik.shape[1])

    total_log_lik = (snp_post * np.log(lik)).sum()
    
    epsilon_prior = beta.logpdf(np.exp(log_epsilon),1,1000)
    skew_prior = beta.logpdf(skew,100,100)

    return(-1*total_log_lik.sum() - epsilon_prior - skew_prior)










