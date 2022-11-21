import autograd.numpy as np
from approxcv.models.abstract_lgiidregr import AbstractiidLGR
from approxcv.acvij import IJ
from approxcv.acvns import NS
from approxcv.models.examples import SLC
import regression.helpers as hp
import os
import time
from pandas import read_csv

args = hp.get_args()
method = "NS"
integrate = False if method == 'IJ' else args.integrate
is_ij = method == "IJ"

modelfolder = "../data-raw/pretrained/slc"
savefolder = "./results/slc"

dat = hp.load_data(modelfolder)
phi = read_csv("%s/phi.csv" % modelfolder).to_numpy()[:,1]
y = hp.approx_y_poisson(np.log(1/phi), dat['y'].to_numpy())
cvfolds = np.arange(dat.shape[0])
J = np.max(cvfolds)

W = read_csv("%s/W.csv" % modelfolder).to_numpy()

prefix = ''
# draws = hp.load_samps(modelfolder, prefix).to_numpy()
(parnames, params) = hp.load_params(modelfolder, prefix)
X = hp.load_x(modelfolder, prefix)[:,0:2]
N,P = X.shape
X = np.concatenate((X, np.diag(np.ones([N]))), axis = 1)
prior = hp.prior_info(modelfolder, prefix)

pm = hp.ParamManager(parnames, params, X).transform_sig()
if integrate:
    reff_mask = ['phi' in x for x in pm.parnames][0:P]
    pm = pm.remove_coeffs(reff_mask).save_params()

offset = dat.log_offset.to_numpy()

model = SLC(pm.X, y, W, prior, pm.parnames, integrate = integrate, log_offset = offset)
acv = NS(model, pm.params, N)

weights = np.ones([N])
model.weighted_loss(pm.params, weights)

pdict = model.unpack_params(pm.params)

prefix_int = "-c" if integrate else ""
with open(savefolder + '/%s_yhats.csv' % (method + prefix_int), 'w') as yhatfile:
    yhatfile.write("i,yhat,loop,loop_time\n")
    for j in range(J + 1):
        print(j)
        o = [i for i in range(N) if cvfolds[i] == j]
        # check for any missing parameters
        if not is_ij:
            acv, model = pm.remove_missing(o).transfer_to_model(acv, model)
            pm = pm.reset_params()
        yhats = hp.compute_yhats(acv, model, o, optim = method == "MCV")
        timedif = yhats['time']
        timedif += ij_time if j == 0 and is_ij else 0
        yhatfile.write(hp.yhat_line('', yhats['yhat'], o, j + 1, timedif))
