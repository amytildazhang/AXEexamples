"""
Latent  Gaussian regression with iid random effects
"""
import autograd.numpy as np
from autograd import grad
import abc
import regression.helpers as hp
ABC = abc.ABC


class AbstractLGR(ABC):
    """ Base class. New models can be defined by inheriting from this class
    """

    def __init__(self, X, Y, prior_info, obs_density, par_names, **kwargs):
        self.N, self.P = X.shape
        self.X = X
        self.Y = Y
        self.R = len(hp.as_list(prior_info['prior_fixef']['location'])) 
        self.prior = prior_info
        self.obs_dist = obs_density
        self.par_names = par_names
        """ Initialize a BuiltinUQ object.
        """

    def weighted_loss(self, params, weights):
        """
        For LOOCV / IF computation within a single sequence. Uses weighted alpha recursion
        :param params:
        :param weights:
        :return:
        """
        param_dict = self.unpack_params(params)
        logp = self.log_prior(param_dict) # log prior: covariance parameters  and fixed effects  
        ll = self.marg_log_obs_lik(param_dict, weights) # scalar, log sum of weighted likelihood of observation model       
        return -logp - ll

    def marg_log_obs_lik(self, pdict, weights):
        phiinv = np.einsum('i,i->i', weights, self.phi_inv(pdict))
        prec = self.marg_prec(pdict, weights)
        mu = self.marg_expect(prec, phiinv, weights)
        ll = hp.ll_normal(np.einsum('i,i->i', self.Y, weights), mu, prec)
        return ll

    def V(self, pdict, weights):
        siginv = self.sigma_inv(pdict)
        phiinv = self.phi_inv(pdict)
        XtPhiinvX = np.einsum("j,ji,j,jk,j->ik", weights, self.X, phiinv, self.X, weights) 
        return np.linalg.inv(XtPhiinvX + siginv)

    def marg_expect(self, prec, phi_inv, weights):
        mu = np.einsum('ij,j,j->i', prec, phi_inv, self.Y)
        return mu 

    def marg_prec(self, pdict, weights):
        V = self.V(pdict, weights)
        X_wt = np.einsum("i,ij->ij", weights, self.X)
        prec = np.einsum("ij,jk,lk->il",  X_wt, V, X_wt)
        return prec

    def predict(self, params, weights):
        pdict = self.unpack_params(params)
        phiinv = self.phi_inv(pdict)
        V = self.V(pdict, weights)
        prec = np.einsum("ij,jk,lk->il",  self.X, V, self.X)
        return self.marg_expect(prec, phiinv, weights)

    def log_prior(self, pdict):
        ll = 0
        for par in pdict.keys():
            if par == 'random': 
                continue
            prior_key = 'prior_' + par
            add_ll = hp.parse_ll(self.prior[prior_key], pdict[par])
            ll = ll + np.sum(add_ll)
        return ll
    
    def unpack_params(self, params):
        masks = self.param_mask()
        param_dict={}
        for par in masks.keys():
            param_dict[par] = np.array([params[i] for i, x in enumerate(masks[par]) if x])
        return param_dict

    # @abc.abstractmethod
    # def log_reff_lik(self, pdict):
    #     raise NotImplementedError

    @abc.abstractmethod
    def sigma_inv(self, *argv, **kwargs):
        # return numpy 1d array 
        raise NotImplementedError

    @abc.abstractmethod
    def phi_inv(self, *argv, **kwargs):
        # return numpy 1d array 
        raise NotImplementedError

    @abc.abstractmethod
    def param_mask(self, *argv, **kwargs):
        # must have keys whicih match keys in prior_dict
        raise NotImplementedError

    


# class AbstractiidLGR(AbstractLGR):
#     def log_reff_lik(self, pdict, weights):
#         """
#         :param beta: 
#         :param sig_params:
#         """
#         if self.cvtype == "LWCV":
#             weights = np.ones([self.N])
#         SigInv = self.sigma_inv(pdict)
#         lrl = hp.ll_normal_ind(pdict['random'], 0, SigInv)
#         wts = self.lrl_weights(weights)
#         weighted_lrl = np.einsum('i,i', lrl, wts)
#         return weighted_lrl




# class AbstractLGmrfR(AbstractLGR):
#     def log_reff_lik(self, pdict, weights):
#         """
#         :param beta: 
#         :param sig_params:
#         """
#         if self.cvtype == "LWCV":
#             weights = np.ones([self.N])
#         SigInv = self.sigma_inv(pdict)
#         lrl = hp.ll_normal(pdict['random'], 0, SigInv)
#         return lrl

# weights = np.ones_like(np.ones([N]))

# pdict = grr.unpack_params(params)
# phiinv = grr.phi_inv(pdict)
# Vmat = V(grr, pdict, weights)
# prec = np.einsum("ij,jk,lk->il",  grr.X, Vmat, grr.X)

# siginv = grr.sigma_inv(pdict)
# phiinv = grr.phi_inv(pdict)
# Vinv = np.einsum("ji,j,jk->ik", grr.X, phiinv, grr.X) + siginv

# np.einsum("ij,i->ij", grr.X, weights)
# Vinv = np.einsum("j,ji,j,jk,j->ik", weights, grr.X, phiinv, grr.X, weights) + siginv


# marg_expect(self, prec, phi_inv, weights):


# def V(grr, pdict, weights):
#     siginv = grr.sigma_inv(pdict)
#     phiinv = grr.phi_inv(pdict)
#     XtPhiinvX = np.einsum("j,ji,j,jk,j->ik", weights, grr.X, phiinv, grr.X, weights) 
#     return np.linalg.inv(XtPhiinvX + siginv)


# phiinv = np.einsum('i,i->i', weights, grr.phi_inv(pdict))
# prec = grr.marg_prec(pdict, weights)
# mu = grr.marg_expect(prec, phiinv, weights)
# ll = hp.ll_normal(np.einsum('i,i->i', grr.Y, weights), mu, prec)
# return ll
