import numpy as np
import pandas as pd
from scipy.stats import beta, binom, betabinom
from scipy.optimize import minimize

def run_snp_mixture_model(nref, nalt, alpha=None, nit=50):
    """ Fit SNP mixture model

    nref -- number of reference bases
    nalt -- number of alternate bases
    alpha -- dirichlet prior on het/homref/homalt/other probabilities
    nit -- number of iterations to run EM algorithm
    """

    if alpha is None:
        alpha = np.array([50, 50, 100, .1])

    # Initialize parameters to be optimized
    pi = alpha / alpha.sum()
    skew = .5
    epsilon = 1e-3

    # Perform EM
    for i in range(0, nit):
        # E-step
        snp_post = method_mixture_postprob(nref, nalt, epsilon=epsilon, skew=skew, pi=pi)

        # M-step
        pi = snp_post.sum(axis=1) + alpha
        pi = pi / pi.sum()

        res = minimize(lambda x: method_mixture_Eloss(snp_post, nref, nalt, x[0], x[1], pi),
                       x0=(np.log(epsilon), skew), bounds=[(np.log(1e-20), np.log(1e-2)), (.4, .6)])

        epsilon = np.exp(res.x[0])
        skew = res.x[1]

        print(f'skew : {skew}, epsilon : {epsilon}, pi : {pi}')

    p_df = pd.DataFrame(snp_post,columns = ['prob_het','prob_homref','prob_homalt','prob_other'])

    return ({'snp_prob':p_df,'epsilon': epsilon, 'skew': skew, 'pi': pi})

# Returns a 4 x N matrix with the posterior probability each SNP is a het, hom_ref, hom_alt, or none of the above
def method_mixture_postprob(nref,nalt,epsilon,skew,pi):

    lik = method_mixture_likelihood(nref,nalt,epsilon=epsilon,skew=skew)
    
    pi_shape = np.ones(len(lik.shape)).astype(int)
    pi_shape[0]=4
    pi = pi.reshape(pi_shape)
    
    post_prob = (lik*pi) / (lik * pi).sum(axis=0)

    return(post_prob)

# Returns a 4 x N matrix with the model likelihood each SNP is a het, hom_ref, hom_alt, or none of the above
def method_mixture_likelihood(nref,nalt,epsilon,skew):
    
    N = nref+nalt
    
    het_lik = binom.pmf(nalt,N,skew)
    homref_lik = binom.pmf(nalt,N,epsilon/3)
    homalt_lik = binom.pmf(nalt,N,1-epsilon/3)
    unif_lik = betabinom.pmf(nalt,N,1,1)
    
    lik = np.stack([het_lik,homref_lik,homalt_lik,unif_lik])
    
    return(lik)
    
# Negative expected log posterior to optimize for the M-step
def method_mixture_Eloss(snp_post,nref,nalt,log_epsilon,skew,pi,min_lik=1e-50):

    lik = method_mixture_likelihood(nref,nalt,np.exp(log_epsilon),skew,pi)
    lik = np.maximum(min_lik,lik)

    total_log_lik = (snp_post * np.log(lik)).sum()
    
    epsilon_prior = beta.logpdf(np.exp(log_epsilon),1,1000)
    skew_prior = beta.logpdf(skew,100,100)

    return(-1*total_log_lik.sum() - epsilon_prior - skew_prior)










