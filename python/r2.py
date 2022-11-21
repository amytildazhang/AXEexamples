from autograd import jacobian
import autograd.numpy as np
from approxcv.models.abstract_lgiidregr import AbstractiidLGR
from approxcv.acvij import IJ
from approxcv.acvns import NS
from approxcv.models.examples import Radon
import regression.helpers as hp
import os
import time
import itertools
from pandas import read_csv

def r2_prefix(mod_no, perc, i, n):
    return "perc%s-n%s-iter%s-mod%s_" % (int(perc*10), n, i, mod_no)

args = hp.get_args()
method = args.method
integrate = False if method == 'MCV' else args.integrate
is_ij = method == "IJ"


modelfolder = "../data-raw/pretrained/r2"
savefolder = "results/r2"

models = [1, 2, 3]
data_perc = [0.3, 0.4, 0.5, 0.6, 0.7]
nc = [3, 4, 6, 9, 12]


loop = 0
combns = itertools.product(models, data_perc, nc, list(range(1, 61)))
prefix_int = "-c" if integrate else ""
with open(savefolder + '/%s_yhats.csv' % (method + prefix_int), 'w') as yhatfile:
    yhatfile.write("model,perc,n_clusters,iter,i,yhat,loop,loop_time\n")
    for x in combns:
        (mod, dperc, n, si) = x
        if si > 35 and dperc == 0.3 and n == 3: 
            continue
        print(x)
        prefix = r2_prefix(mod, dperc, si, n)
        (parnames, params) = hp.load_params(modelfolder, prefix)
        dat = read_csv("%s/%sdata.csv" % (modelfolder, prefix))
        y = dat['log_radon'].to_numpy()
        X = hp.load_x(modelfolder, prefix)
        N,P = X.shape
        X[:,0] = 1.0
        prior = hp.prior_info(modelfolder, prefix)
        # pm_full = hp.optim_fulldata(parnames, params, X, Radon, y, prior)

        pm = hp.ParamManager(parnames, params, X)
        # pm = hp.ParamManager(pm_full.parnames, pm_full.params, pm_full.X)
        if integrate:
            county_mask = ['county' in x for x in pm.parnames][0:P]
            pm = pm.remove_coeffs(county_mask).save_params()

        model = Radon(pm.X, y, prior, pm.parnames,  integrate = integrate)
        ij_start = time.time()
        acv = IJ(model, pm.params, N) if is_ij else NS(model, pm.params, N)
        ij_time = time.time() - ij_start

        o = [i for (i, x) in enumerate(dat['county'].to_list()) if x == 'OLMSTED']

        if not is_ij:
            acv, model = pm.remove_missing(o).transfer_to_model(acv, model)
        yhats = hp.compute_yhats(acv, model, o, optim = method == "MCV")
        pm = pm.reset_params()

        timedif = yhats['time']
        timedif += ij_time if is_ij else 0

        y_prefix = "%d,%1f,%d,%d," % (mod, dperc, n, si)
        yhatfile.write(hp.yhat_line(y_prefix, yhats['yhat'], np.arange(23) + 1, loop, timedif))
        loop += 1
