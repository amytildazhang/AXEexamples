"""
ACV with Newton-Raphson step
"""
import autograd
import autograd.numpy as np
from scipy.optimize import minimize


class NS:
    def __init__(self, model, theta_one, T):
        """
        :param model: Must implement weighted_loss
        :param theta_one: MLE / MAP fit on the entire dataset (with all weights set to 1)
        :param T: length of sequence / number of sites in a MRF
        :param o: Indices of sites in the held out fold
        """
        self.weights_one = np.ones([T])
        self.theta_one = theta_one
        self.loss_fun = model.weighted_loss
        # self.H = self.compute_hessian(self.weights_one)
    
    def optimize(self, o = [], save = False):
        weights = np.ones_like(self.weights_one)
        weights[o] = 0  # dropped sites
        params = minimize(self.loss_fun, 
            self.theta_one, method='Newton-CG',
            jac = self.compute_Jvec, hess = self.compute_hessian,
            args = (weights))
        if save:
            self.theta_one = params['x']
        return params['x']
    
    def compute_hessian(self, params, weights):
        eval_hess = autograd.hessian(self.loss_fun, argnum=0)
        H = eval_hess(params, weights)
        return H
    
    def compute_Jvec(self, params, weights):
        eval_jvec = autograd.jacobian(self.loss_fun, argnum=0)
        J = eval_jvec(params, weights)
        return J
    
    def compute_params_acv(self, o):
        weights = np.ones_like(self.weights_one)
        weights[o] = 0  # dropped sites
        Hinv = np.linalg.inv(self.compute_hessian(self.theta_one, weights))
        # HdeltaInv = np.linalg.inv(self.H - self.compute_hessian(1 - weights))
        params_acv = self.theta_one + Hinv.dot(self.compute_Jvec(self.theta_one, weights))
        return params_acv
