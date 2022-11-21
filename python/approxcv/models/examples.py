from approxcv.models.abstract_lgiidregr import AbstractiidLGR, PoissonLGMRF
import autograd.numpy as np
import regression.helpers as hp
from autograd.scipy.linalg import block_diag

class Radon(AbstractiidLGR):
    """
    """
    def __init__(self, X, Y, prior_info, par_names, integrate):
        # obs_density = {'dist': 'normal', 'invlink': lambda x: x}
        super().__init__(X=X, Y=Y, prior_info=prior_info, 
         par_names=par_names,  integrate = integrate)
        """
        """

    def sigma(self, pdict):
        sig = np.array(pdict['covariance'])
        if not self.integrate:
            P2 = len(pdict['random'])
            sig = np.ones([P2])*sig
        return sig
    
    def phi(self, pdict):
        phi_val = pdict['aux']**2
        if self.integrate:
            phi_val = phi_val + self.sigma(pdict)
        return np.ones([self.N])*np.array(phi_val)
    
    def param_mask(self):
        parnames = self.par_names
        aux_mask = ['sigma' in x for x in parnames]
        random_mask = ['county' in x and 'Sigma' not in x for x in parnames]
        covariance_mask = ['Sigma' in x for x in parnames]
        fixef_mask = [not x and not y and not z for (x,y,z) in zip(aux_mask, covariance_mask,random_mask)]
        pdict = {'aux': aux_mask, 'covariance':covariance_mask, 'fixef':fixef_mask}
        if not self.integrate:
            pdict['random'] = random_mask
        return pdict




class Eight(AbstractiidLGR):
    """
    """
    def __init__(self, X, Y, prior_info, par_names, integrate):
        super().__init__(X=X, Y=Y, 
            prior_info=prior_info,  par_names=par_names, integrate = integrate)
        """
        """
    def log_prior(self, pdict):
        if self.integrate:
            ll = 0
        else:
            ll = self.log_reff(pdict)
        return ll

    def sigma(self, pdict):
        sig = np.array(pdict['covariance'])**2
        P2 = pdict['random'].shape[0]
        return  np.ones([P2])*sig
            
    def phi(self, pdict):
        phi_val = np.array([15.0, 10.0, 16, 11, 9, 11, 10, 18])**2
        if self.integrate:
            phi_val = phi_val + self.sigma(pdict)[0]
        return phi_val
    
    def param_mask(self):
        parnames = self.par_names
        random_mask = ['theta' in x for x in parnames]
        covariance_mask = ['tau' in x for x in parnames]
        fixef_mask = [x == 'mu' for x in parnames]
        pdict = {'covariance':covariance_mask, 'fixef':fixef_mask}
        if not self.integrate:
            pdict['random'] = random_mask
        return pdict



class LOL(AbstractiidLGR):
    """
    """
    def __init__(self, X, Y, prior_info, par_names, integrate):
        super().__init__(X=X, Y=Y, 
            prior_info=prior_info,  par_names=par_names, integrate = integrate)
        """
        """
    def predict(self, params):
        mu = super().predict(params)
        return np.exp(mu)

    def log_obs(self, pdict):
        if self.integrate:
            ll = super().log_obs(pdict)
        else:
            mu = self.expectation(pdict)
            ll = hp.ll_poisson(self.Y, np.exp(mu), 1.0)
        return ll

    def sigma(self, pdict):
        n_champ = np.sum(['champion' in x for x in self.par_names]) - 1
        sigs = pdict['covariance'][1]*np.ones([n_champ])
        if not self.integrate:
            n_player = np.sum(['player' in x for x in self.par_names]) - 1
            sigs_player = pdict['covariance'][1]*np.ones([n_player])
            sigs = np.concatenate((sigs, sigs_player))
        return sigs
    
    def phi(self, pdict):
        phis = np.exp(self.expectation(pdict))
        if self.integrate:
            phis = phis + np.array(pdict['covariance'][1])
        return np.array(phis)
    
    def param_mask(self):
        pc_mask = ['champion' in x or 'player' in x for x in self.par_names]
        covariance_mask = ['Sigma' in x for x in self.par_names]
        random_mask = [x and not y for (x,y) in zip(pc_mask, covariance_mask)]
        fixef_mask = [not x for x in pc_mask]
        self.R = sum(random_mask) + sum(fixef_mask)
        return {'random': random_mask, 'covariance':covariance_mask, 'fixef':fixef_mask}


class SLC(PoissonLGMRF):
    """
    """
    def __init__(self, X, Y, W, prior_info, par_names, integrate, log_offset, 
        remove_coeff=None):
        # obs_density = {'dist': 'normal', 'invlink': lambda x: x}
        super().__init__(X=X, Y=Y,W=W, prior_info=prior_info, 
         par_names=par_names, integrate = integrate, log_offset = log_offset)
        """
        """
    # def log_reff(self, pdict):
    #     reffs = pdict['random']
    #     if not self.o is None:
    #         reffs = np.insert(reffs, self.o, 0)
    #     siginv = self.sigma(pdict)
    #     ll_reffs = hp.ll_normal(pdict['random'], 0.0, siginv)
    #     return ll_reffs
    
    def sigmainv(self, pdict): 
        D = np.diag(self.D)
        siginv = pdict['tau'] * (D - pdict['alpha'] * self.W)
        return siginv

    def param_mask(self):
        parnames = self.par_names
        fixef_mask = ['beta' in x for x in parnames]
        alpha_mask = [x == 'alpha' for x in parnames]
        tau_mask = [x == 'tau' for x in parnames]
        pdict = {'alpha': alpha_mask, 'tau':tau_mask, 'fixef':fixef_mask}
        if not self.integrate:
            pdict['random'] = ['phi' in x for x in parnames]
        return pdict

    # def predict(self, params, o = None):
    #     mu = super().predict(params, o)
    #     if not self.integrate and o is not None:
    #         mu = np.log(mu - self.log_offset)
    #         weights = np.ones(self.N)
    #         weights[o] = 0
    #         nbhd = self.W[o, :][:, weights == 1]
    #         pdict = self.unpack_params(params)
    #         alpha = np.einsum('ij,j', nbhd, pdict['random'])/np.sum(nbhd)
    #         mu = np.exp(mu + alpha + self.log_offset)
    #     return mu







class AIR(PoissonLGMRF):
    """
    """
    def __init__(self, X, Y, W, prior_info, par_names, integrate, log_offset):
        # obs_density = {'dist': 'normal', 'invlink': lambda x: x}
        # obs_density = {'dist': 'normal', 'invlink': lambda x: x}
        super().__init__(X=X, Y=Y,W=W, prior_info=prior_info, 
         par_names=par_names, integrate = integrate, log_offset = log_offset)
        self.J = W.shape[0] # number of IGs
        self.T = int(X.shape[0]/self.J) # timepoints
        """
        """
    def sigmainv(self, pdict):
        H = self.H_sigma()
        ones = np.ones([self.N])
        IrhH = np.diag(np.ones_like(ones)) - pdict['rho_t'] * H
        Sigma = self.Q_sigma(pdict['rho_s'])
        HSigma = np.dot(IrhH, Sigma)
        tauinv = pdict['tau2']**-1
        return np.einsum('ij,jk->ik', HSigma, IrhH)*tauinv

    def param_mask(self):
        parnames = self.par_names
        rho_s = ['rho_s' == x for x in parnames]
        rho_t = ['rho_t' == x for x in parnames]
        tau = ['tau2' == x for x in parnames]
        fixef_mask = ['beta' in x for x in parnames]
        pdict = {'fixef':fixef_mask, 'rho_s': rho_s, 'rho_t': rho_t, 'tau2': tau}
        if not self.integrate:
            pdict['random'] = ['phi' in x for x in parnames]
        return pdict
    
    def H_sigma(self):
        P = self.J * self.T
        ones = np.ones([self.J * (self.T - 1)])
        H = np.zeros((P, P))
        H[self.J:P, :][:, 0:(P - self.J)] = np.diag(ones)
        return H
    
    def block_row(self, Sigma, t):
        zeros_l = np.zeros((self.J, self.J * t))
        zeros_r = np.zeros((self.J, self.J * (self.T - t - 1)))
        return np.concatenate((zeros_l, Sigma, zeros_r), axis = 1)
    
    def block_diag(self, Sigma):
        rows = (self.block_row(Sigma, t) for t in np.arange(self.T))
        return np.concatenate(rows, axis = 0)
    
    def Q_sigma(self, rhoj):
        ones = np.ones([self.J])
        Q = rhoj * (np.diag(self.D) - self.W) + (1 - rhoj)*np.diag(ones)
        Sigma = Q
        return self.block_diag(Sigma)

# mat = np.concatenate((block_row(model, Sigma, t) for t in np.arange(model.T)), axis = 0)

# def block_row(model, Sigma, t):
#     row = Sigma
#     zeros_l = np.zeros((model.J, model.J * t))
#     zeros_r = np.zeros((model.J, model.J * (model.T - t - 1)))
#     return np.concatenate((zeros_l, Sigma, zeros_r), axis = 1)


# P = model.J * model.T
# ones = np.ones([model.J * (model.T - 1) + 1])
# H = np.zeros((P, P))
# H[J:P, :][:, 0:(P - J)] = np.diag(ones)
 

# pdict = model.unpack_params(pm.params)
# mu = model.expectation(pdict)
# Phi = model.phi(pdict)
# weights = np.ones([model.N])
# weights[o] = 0
# train_mask = weights == 1
# test_mask = weights != 1
# Sigtrain_inv = np.linalg.inv(Phi[train_mask, :][:, train_mask])
# Sigttrain = Phi[test_mask, :][:, train_mask]
# dif = np.einsum('ij,jk,k->i', Sigttrain, Sigtrain_inv, (model.Y[train_mask] - mu[train_mask]))
# mu_test = np.exp(mu[test_mask] + dif + model.log_offset[test_mask])

# def log_reff(self, pdict):
#     # reffs = pdict['random']
#     # rhoj = pdict['']
#     # rhot = pdict['']
#     # sigma = pdict['']
#     # siginv = self.Q_sigma(pdict)
#     # ll_reffs = 0
#     # if not self.o is None:
#     #     reffs = np.insert(reffs, self.o, 0)
#     # for time in np.arange(self.T):
#     #     reff_t = self.reffs_at_t(time, reffs)
#     #     mu = 0 if time == 0 else reffs_at_t(time - 1, reffs)
#     #     ll_reffs = ll_reffs + hp.ll_normal(reff_t, mu, siginv * sigma**-2)
#     # siginv = self.sigma(pdict)
#     # ll_reffs = hp.ll_normal(pdict['random'], 0.0, siginv)
#     return 0
#     # def reffs_at_t(time, reffs):
#     #     return np.array(reffs[time*self.J:(time + 1)*self.J])
