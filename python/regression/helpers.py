import autograd.numpy as np
import scipy.stats as stats
from pandas import read_csv
import json
from copy import deepcopy
import time
import argparse
from approxcv.acvns import NS

# for Poisson GLMMs, use working response as in Breslow and Clayton
def approx_y_poisson(eta, y):
    mu = np.exp(eta)
    return eta + (y - mu) / mu

def posterior_mode(vec):
    kernel = stats.gaussian_kde(vec)
    height = kernel.pdf(vec)
    mode_value = vec[np.argmax(height)]
    return mode_value

def param_mode(matrix):
    modes = [posterior_mode(matrix[:,i]) for i in range(matrix.shape[1])]
    return np.array(modes)

def posterior_mean(vec):
    return np.mean(vec)

def param_mean(matrix):
    modes = [posterior_mean(matrix[:,i]) for i in range(matrix.shape[1])]
    return np.array(modes)

# functions for getting logliks 

def lldict_scale(ll_dict):
    if ll_dict['adjusted_scale'] is None:
        scale = ll_dict['scale']
    else:
        scale = ll_dict['adjusted_scale']
    return scale


def pack_ll(dist, par1, par2 = None):
    lldict = {'dist': dist}
    if dist == 'normal':
        lldict['location'] = par1
        lldict['scale'] = par2
        lldict['adjusted_scale'] = None
    elif dist == 'exponential':
        lldict['scale'] = par1
        lldict['adjusted_scale'] = None
    elif dist == 'decov':
        lldict['shape'] = par1
        lldict['scale'] = par2
    elif dist == 'poisson':
        lldict['rate'] = par1
        lldict['offset'] = par2 if par2 is not None else 1
    return lldict

def parse_ll(ll_dict, x):
    dist = ll_dict['dist']
    if dist == 'normal':
        std = np.array(lldict_scale(ll_dict))*1.0
        ll = ll_normal_ind(x, np.array(ll_dict['location']), std**-2)
    elif dist == 'exponential':
        ll = ll_exponential(x, np.array(lldict_scale(ll_dict)))
    elif dist == 'decov' or dist == 'gamma':
        ll = ll_gamma(x, np.array(ll_dict['shape']), np.array(ll_dict['scale']))
    elif dist == 'poisson':
        ll = ll_poisson(x, np.array(ll_dict['rate']), np.array(ll_dict['offset']))
    elif dist == 'invgamma':
        ll = ll_invgamma(x, np.array(ll_dict['shape']), np.array(ll_dict['scale']))
    elif dist == 'uniform':
        ll = 0
    return ll

def ll_gamma(x, shape, scale):
    ll = -shape * np.log(scale) + (shape - 1)*np.log(x) - x/scale
    return np.array(ll)

def ll_invgamma(x, shape, scale):
    ll = shape * np.log(scale) - (shape - 1)*np.log(x) - scale/x
    return np.array(ll)


def ll_exponential(x, scale):
    ll = -np.log(scale) - x/scale
    return np.array(ll)
    # ll_exponential(1, 1/3) should be -1.9

def ll_poisson(x, rate, offset):
    """
    x: numpy 1darray
    rate: numpy 1darray
    """
    rate = rate * offset
    log_factorial = np.array([np.sum(np.log(np.arange(i) + 1)) for i in x])
    ll = x * np.log(rate) - rate - log_factorial
    return np.array(ll)

def ll_normal_ind(x, mu, siginv): # assumes siginv is a vector of same dimensions as x
    ll = -0.5 * siginv *( (x -  mu))**2 + 0.5 * np.log(siginv)     
    return np.array(ll)

def ll_normal(x, mu, siginv): # assumes siginv is positive-definite
    sign, logdet = np.linalg.slogdet(siginv)
    ll = -0.5 * np.einsum('i,ij,j', x - mu, siginv, x - mu) + 0.5 * sign * logdet
    return np.array(ll)

def as_list(val):
    return val if isinstance(val, list) else [val]

def merge_prior(prior_dict):
    fixef = {}
    prior = prior_dict['prior']
    intercept = prior_dict['prior_intercept']
    if not prior_dict['prior'] is None:
        fixef['location'] = as_list(intercept['location']) + as_list(prior['location'])
        fixef['adjusted_scale'] = as_list(lldict_scale(intercept)) + as_list(lldict_scale(prior))
        fixef['scale'] = None
        fixef['dist'] = 'normal'
    else:
        fixef =  prior_dict['prior_intercept']
    prior_dict.pop('prior')
    prior_dict.pop('prior_intercept')
    prior_dict['prior_fixef'] = fixef
    return prior_dict


# functions for reading in posterior info from pretrained models


def acv_yhats(acv, grr, o):
    check_missing_X(X, o, )

def prior_info(dirpath, prefix=''):
    prior_path = "%s/%sprior.json" % (dirpath, prefix)
    with open(prior_path, 'r') as prior_file:
        prior = json.load(prior_file)
    return merge_prior(prior)

def load_samps(dirpath, prefix):
    fpath = "%s/%sallsamps.csv" % (dirpath, prefix)
    df = read_csv(fpath)
    return df

def load_params(dirpath, prefix):
    fpath = "%s/%smap.csv" % (dirpath, prefix)
    df = read_csv(fpath)
    params = df.map.to_numpy().flatten()
    parnames = df.param.to_list()
    return (parnames, params)

def load_x(dirpath, prefix):
    fpath = "%s/%sX.csv" % (dirpath, prefix)
    X = read_csv(fpath).to_numpy()
    return X

def load_w(dirpath, prefix):
    fpath = "%s/%sW.csv" % (dirpath, prefix)
    X = read_csv(fpath).to_numpy()
    return X

def loglik(params, model):
    pdict = model.unpack_params(params)
    return model.marg_log_obs_lik(pdict)

def load_data(dirpath, prefix=''):
    fpath = "%s/%sdata.csv" % (dirpath, prefix)
    dat = read_csv(fpath)
    return dat

# def param_names(dirpath, prefix):
#     df = load_samps(dirpath, prefix)
#     return df.columns.tolist()

class ParamManager:
    def __init__(self, parnames, params, X):
        self.X = deepcopy(X)
        self.params = deepcopy(params)
        self.parnames = deepcopy(parnames)
        self.og = deepcopy((parnames, params, X))
        self.sigvars = ['Sigma', 'sigma', 'tau', 'tau2', 'alpha', 'rho_t', 'rho_s']
     
    def coeff_idx(self):
        coeff_mask = [not any([y in x for y in ['Sigma', 'sigma', 'tau']]) for x in self.parnames]
        return self.mask_to_idx(coeff_mask)
    
    def return_params(self):
        return (self.parnames, self.params, self.X)
    
    @staticmethod
    def mask_to_idx(mask):
        return [i for (i,x) in enumerate(mask) if x]
    
    def remove_coeffs(self, exclude_mask):
        sig_idx = self.mask_to_idx([any(y in x for y in self.sigvars) for x in self.parnames])
        keep_coeff_idx = self.mask_to_idx([not x for x in exclude_mask]) 
        keep_params = keep_coeff_idx + sig_idx
        self.params = self.params[keep_params]
        self.parnames = [self.parnames[i] for i in keep_params]
        self.X = self.X[:, keep_coeff_idx]
        return self
    
    def remove_missing(self, o):
        N = self.X.shape[0]
        train_idx = [i for i in range(N) if i not in o]
        colsums = np.einsum("ij->j", self.X[train_idx, :])
        return self.remove_coeffs([x == 0 for x in colsums])
    
    def transform_sig(self):
        sig_mask = [any([y == x for y in self.sigvars]) for x in self.parnames]
        sig_params = np.log(self.params[sig_mask])
        self.params[sig_mask] = sig_params
        return self
    
    def transfer_to_model(self, acv, model):
        model.X = self.X
        model.par_names = self.parnames
        acv.theta_one = self.params
        return (acv, model)
    
    def reset_params(self):
        self.parnames, self.params, self.X = deepcopy(self.og)
        # self.parnames = parnames
        # self.params = params
        # self.X = X
        return self

    def save_params(self):
        self.og = deepcopy((self.parnames, self.params, self.X))
        return self

def compute_yhats(acv, model, o, optim = False):
    start = time.time()
    if optim:
        params = acv.optimize()
    else:
        params = acv.compute_params_acv(o)
    yhats = model.predict(params)[o]
    timedif = time.time() - start
    return {'yhat': yhats, 'time': timedif}

def optim_fulldata(parnames, params, X, modelfun, y, prior, **kwargs):
    pm = ParamManager(parnames, params, X).transform_sig()
    model = modelfun(X = pm.X, Y = y, prior_info = prior, 
        par_names = pm.parnames, integrate = False, **kwargs)
    acv = NS(model, pm.params, pm.X.shape[0])
    params_map = acv.optimize()
    pm.params = params_map
    return pm.save_params()


def yhat_line(prefix, yhats, o, loop, time):
    o_list = [e + 1 for e in o]
    all_str = zip(o_list, yhats)
    yhatlines = [','.join([str(y) for y in x]) for x in all_str]
    looptime = ','.join([str(loop),str(time)])
    yhatlines = [prefix + x + "," + looptime for x in yhatlines]
    return '\n'.join(yhatlines) + '\n'


# def param_line(prefix, parnames, params, params_acv, loop):
#     param_str = [str(e) for e in params]
#     acv_str = [str(e) for e in params_acv]
#     lines = [prefix + ','.join([x,y,z,str(loop)]) for x,y,z in zip(parnames, param_str, acv_str)]
#     return '\n'.join(lines) + '\n'

def get_args():
    parser = argparse.ArgumentParser(description='Set ACV method.')
    parser.add_argument('method', help='ACV method. Choose one of IJ, NS, or MCV.', choices=['IJ', 'NS', 'MCV'])
    parser.add_argument('--integrate',  help = 'Whether to integrate out theta_j',
        action = 'store_true')
    args = parser.parse_args()
    return args

# sig_mask = [any([y in x] for y in pm.sigvars) for x in pm.parnames]
# sig_params = np.log(pm.params[sig_mask])
# pm.params[sig_mask] = sig_params

# sig_idx = pm.mask_to_idx(['sigma' in x or 'Sigma' in x for x in pm.parnames])
# keep_coeff_idx = pm.mask_to_idx([not x for x in exclude_mask]) 
# keep_params = keep_coeff_idx + sig_idx
# pm.params = pm.params[keep_params]
# pm.parnames = [pm.parnames[i] for i in keep_params]
# pm.X = pm.X[:, keep_coeff_idx]
# # return self

