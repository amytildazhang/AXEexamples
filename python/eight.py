import autograd.numpy as np
from approxcv.models.abstract_lgiidregr import AbstractiidLGR
from approxcv.acvij import IJ
from approxcv.acvns import NS
from approxcv.models.examples import Eight
import regression.helpers as hp
import os
import time

# possible values for method: NS, IJ, MCV
args = hp.get_args()
method = args.method
integrate = False if method == 'MCV' else args.integrate
is_ij = method == "IJ"

modelfolder = "../data-raw/pretrained/eight_schools"
savefolder = "./results/eight_schools"

J = 8
y = np.array([28.0, 8, -3, 7, -1, 1, 18, 12])
sd = np.array([15.0, 10, 16, 11, 9, 11, 10, 18])
scale_array = np.round(np.arange(0.1, 4.1, step = 0.1), decimals = 1)
ones = np.ones([J])
X = np.hstack((ones[:,np.newaxis], np.diag(ones)))

prior = {'prior_covariance': {'dist': 'uniform'},'prior_fixef':{'dist':'uniform', 'location':[0]}}
cvfolds = range(8)
P = J + 1

prefix_int = "-c" if integrate else ""
with open(savefolder + '/%s_yhats.csv' % (method + prefix_int), 'w') as yhatfile:
    yhatfile.write("data_scale,i,yhat,loop,loop_time\n")
    for scl in scale_array:
        print(scl)
        prefix = '%.1f_' % scl
        (parnames, params) = hp.load_params(modelfolder, prefix)
        ordered_parnms = [x for x in parnames if x != "tau"] + ['tau']
        reorder = [np.where([parnames[i] == x for (i, z) in enumerate(parnames)])[0].tolist()[0] for x in ordered_parnms]
        ordered_params = params[reorder]

        # pm_full = hp.optim_fulldata(parnames, params, X, Eight, y, prior)
        loop_no = 0
        pm = hp.ParamManager(ordered_parnms, ordered_params, X).transform_sig()
        # pm = hp.ParamManager(pm_full.parnames, pm_full.params, pm_full.X)
        if integrate:
            reff_mask = ['theta' in x for x in pm.parnames][0:J]
            pm = pm.remove_coeffs(reff_mask).save_params()

        model = Eight(pm.X, scl*y, prior, pm.parnames,  integrate = integrate)
        ij_start = time.time()
        acv = IJ(model, pm.params, J) if is_ij else NS(model, pm.params, J)
        ij_time = time.time() - ij_start
        # save gradient
        # save hessian   
        for j in cvfolds:
            print(j)
            loop_no += 1
            o = [j]
            # check for any missing parameters
            if not is_ij:
                acv, model = pm.remove_missing(o).transfer_to_model(acv, model)
            yhats = hp.compute_yhats(acv, model, o, optim = method == "MCV")
            pm = pm.reset_params()

            timedif = yhats['time']
            timedif += ij_time if j == 0 and is_ij else 0
            yhatfile.write(hp.yhat_line('%.1f,' % scl, yhats['yhat'], o, loop_no, timedif))
