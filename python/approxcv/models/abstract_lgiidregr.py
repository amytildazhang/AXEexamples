"""
Latent  Gaussian regression with independent Y_j | x_j'beta, phi_j 
"""
import autograd.numpy as np
from autograd import grad
import abc
import regression.helpers as hp
ABC = abc.ABC


#"B" weighting type
class AbstractiidLGR(ABC):
    """ Base class. New models can be defined by inheriting from this class
    """

    def __init__(self, X, Y, prior_info, par_names, integrate, **kwargs):
        self.N, self.P = X.shape
        self.X = X
        self.Y = Y
        self.R = len(hp.as_list(prior_info['prior_fixef']['location'])) # set in unpack_params
        self.prior = prior_info
        # self.obs_dist = obs_density
        self.par_names = par_names
        self.integrate = integrate
        self.sigvars = ['covariance', 'aux',  'tau2', 'tau', 'alpha', 'rho_t', 'rho_s']
        """ Initialize a BuiltinUQ object.
        """

    def predict(self, params):
        param_dict = self.unpack_params(params)
        mu = self.expectation(param_dict)
        return mu

    def weighted_loss(self, params, weights):
        """
        For LOOCV / IF computation within a single sequence. Uses weighted alpha recursion
        :param params:
        :param weights:
        :return:
        """
        param_dict = self.unpack_params(params)
        logp = self.log_prior(param_dict) # log prior: covariance parameters  and fixed effects  
        ll = self.log_obs(param_dict) # scalar, log sum of weighted likelihood of observation model       
        return -logp - np.einsum('i,i', ll, weights)

    def log_reff(self, pdict):
        siginv = self.sigma(pdict)**-1
        ll_reffs = hp.ll_normal_ind(pdict['random'], 0.0, siginv)
        return np.sum(ll_reffs)

    def log_obs(self, pdict):
        prec = self.prec(pdict)
        mu = self.expectation(pdict)
        ll = hp.ll_normal_ind(self.Y, mu, prec)
        return ll

    def expectation(self, pdict):
        if 'random' in pdict.keys():
            effs = np.concatenate((pdict['fixef'], pdict['random']))
        else:
            effs =pdict['fixef']
        # X = self.get_X()
        mu = np.einsum('ij,j->i', self.X, effs)
        return mu

    def prec(self, pdict):
        return self.phi(pdict)**-1 # in the iid case, X2SigmaX2' is sigmaI

    def log_prior(self, pdict):
        ll = 0
        for par in pdict.keys():
            prior_key = 'prior_' + par
            if par == 'random': 
                add_ll = self.log_reff(pdict)
            else:
                add_ll = hp.parse_ll(self.prior[prior_key], pdict[par])
            ll = ll + np.sum(add_ll)
        return ll
    # 
    # def get_X(self):
    #     X = self.X
    #     return X

    def unpack_params(self, params):
        masks = self.param_mask()
        param_dict={}
        for par in masks.keys():
            param_dict[par] = np.array([params[i] for i, x in enumerate(masks[par]) if x])
            if par in self.sigvars:
                param_dict[par] = np.exp(param_dict[par])
        return param_dict

    # @abc.abstractmethod
    # def log_reff_lik(self, pdict):
    #     raise NotImplementedError

    @abc.abstractmethod
    def sigma(self, *argv, **kwargs):
        # return numpy 1d array 
        raise NotImplementedError

    @abc.abstractmethod
    def phi(self, *argv, **kwargs):
        # return numpy 1d array 
        raise NotImplementedError

    @abc.abstractmethod
    def param_mask(self, *argv, **kwargs):
        # must have keys whicih match keys in prior_dict
        raise NotImplementedError

 

class PoissonLGMRF(AbstractiidLGR):
    def __init__(self, X, Y, W, prior_info, par_names, integrate, log_offset = 0.0):
        # obs_density = {'dist': 'normal', 'invlink': lambda x: x}
        super().__init__(X=X, Y=Y, prior_info=prior_info, 
         par_names=par_names,  integrate = integrate)
        self.o = 0
        self.W = W
        self.D = np.einsum('ij->i', W)
        self.log_offset = log_offset
        self.keep_mask = [True for x in Y]

    def filter_matrix(self, matrix):
        return matrix[self.keep_mask, :][:, self.keep_mask]

    def log_reff(self, pdict):
        if self.integrate:
            ll = 0
        else:
            reffs = pdict['random']#np.einsum('i,i->i', self.keep_mask, pdict['random'])
            siginv = self.filter_matrix(self.sigmainv(pdict))
            ll = hp.ll_normal(reffs, 0.0, siginv)
        return ll

    def sigmainv(self, pdict):
        return 0

    def param_mask(self):
        return None

    def phi(self, pdict):
        phi = np.exp(self.expectation(pdict))
        if self.integrate:
            phi = np.diag(phi) + self.sigma(pdict)
        return phi

    def sigma(self, pdict):
        siginv = self.sigmainv(pdict)
        return np.linalg.inv(siginv)

    def prec(self, pdict):
        phi = self.phi(pdict)
        return np.linalg.inv(phi)

    def log_obs(self, pdict):
        mu = self.expectation(pdict)
        if self.integrate:
            prec = self.filter_matrix(self.prec(pdict))
            ll = hp.ll_normal(self.Y[self.keep_mask], mu[self.keep_mask], prec)
        else:
            ll = hp.ll_poisson(self.Y, np.exp(mu), 1.0)
        return ll

    def cond_exp(self, mu, Sigma):
        out_mask = [not x for x in self.keep_mask]
        keep_siginv = np.linalg.inv(self.filter_matrix(Sigma))
        out_to_keep = Sigma[out_mask, :][:, self.keep_mask]
        delta = np.einsum('ij,jk,k->i', out_to_keep, keep_siginv, mu[self.keep_mask])
        return mu[out_mask] +  delta

    def predict(self, params, o = None):
        mu = super().predict(params)
        pdict = self.unpack_params(params)
        weights = np.ones([self.N])
        weights[o] = 0
        if self.integrate:
            self.keep_mask = weights == 1
            mu_test = self.cond_exp(mu, self.phi(pdict))
        else:
            test_mask = weights != 1    
            mu_test = mu[test_mask]
        return np.exp(mu_test + self.log_offset)

    def weighted_loss(self, params, weights):
        """
        For LOOCV / IF computation within a single sequence. Uses weighted alpha recursion
        :param params:
        :param weights:
        :return:
        """
        # if self.integrate:
        self.keep_mask = weights == 1.0
        # else:
        # self.keep_mask = weights
        param_dict = self.unpack_params(params)
        logp = self.log_prior(param_dict) # log prior: covariance parameters  and fixed effects  
        ll = self.log_obs(param_dict) # scalar, log sum of weighted likelihood of observation model       
        if not self.integrate:
            ll = np.einsum('i,i', ll, weights)
        return -logp - ll

