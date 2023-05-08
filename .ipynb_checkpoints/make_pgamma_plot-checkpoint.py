# ratkiewicz1@llnl.gov
import pandas as pd
import uproot
import numpy as np
import os

from datetime import datetime
from matplotlib import pyplot as plt
from matplotlib import colors
from io import StringIO
from lmfit import Minimizer, Parameters, report_fit, fit_report
import scipy
import math

def estimate_parameters(xVals, yVals):
    # does some simple parameter estimation that might be OK.
    # need:
    # sigma
    # volume term
    pars = {}
    pars['A'] = max(yVals)
    pars['sigma'] = 1.2 # just a guess
    # we will pass this two extra bins on either side
    # y1 = m*x1 + b
    # y2 = m*x2 + b
    # --> y1-y2 = m(x1 - x2) --> m = (y1 - y2)/(x1 - x2)
    # b = y1 - x1*(y1-y2)/(x1 - x2)

    #pars['bg_slope']  = (yVals[2] - yVals[-2])/(xVals[2] - xVals[-2])
    #if pars['bg_slope'] == 0.0:
    #    pars['bg_slope'] = 0.01
    #pars['bg_offset'] = yVals[2] - xVals[2]*pars['bg_slope']
    #pars['step_amp']  = (yVals[2]+yVals[1])*0.25*0.5 # guess
    #if pars['step_amp'] == 0.0:
    #   pars['step_amp'] = 0.01
    return pars

def find_nearest(array,value):
    # finds element and index closest to value
    # from: https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) \
            < math.fabs(value - array[idx])):
        return array[idx-1],idx-1
    else:
        return array[idx],idx

def gaussian(params,e_gamma):
    # standard gaussian for fit
    return params['A']*np.exp(-0.5*pow((e_gamma - params['mean'])/params['sigma'],2.))

def linear_bg(params,e_gamma):
    # standard linear bg
    return params['bg_offset']
    #+params['bg_slope']*e_gamma

def step_bg(params,e_gamma):
    #return params['step_amp']/pow(1. + np.exp(e_gamma),2.)
    # y = constant * ERFC( (x-c)/(SQRT(2)*sigma) )
    # this is the radware step function: https://radware.phy.ornl.gov/gf3/gf3.html#1.
    return params['step_amp']*scipy.special.erfc((e_gamma - params['mean'])\
                /(np.sqrt(2.)*params['sigma']))

def fit_gamma(Params,e_gamma, data, dataUnc):
    model = gaussian(Params, e_gamma) +\
            linear_bg(Params, e_gamma) #+\
            #step_bg(Params, e_gamma)

    resids = model - data
    weighted = np.sqrt(pow(resids,2.0)/pow(dataUnc,2.0))

    return weighted

in_data = {}
with uproot.open('output.root') as file:
    in_data['ExEg'] = file['hExc_v_Egam;1'] # get excitation energy matrix
    in_data['Singles'] = file['hExc_Subtracted_Sum;1']   # get particle-singles spectrum

# now there is a numpy matrix with the data in it, so we are free from the 
# tyranny and dumpster fire of memory management that is root.

# this assumes gamma energy on the X axis, which (while not incorrect) is heathenish.
# if you (righteously) have gammas on the Y axis switch the indicies.

exen_index  = 2
gamma_index = 1
data_index = 0
CN = '96Mo'
data = {}

data['ExEg'] = {'matrix' : in_data['ExEg'].to_numpy()}
data['ExEg']['gamma axis'] = data['ExEg']['matrix'][gamma_index]
data['ExEg']['exeg axis'] = data['ExEg']['matrix'][exen_index]

data['ExEg']['gamma centers'] = [(data['ExEg']['matrix'][gamma_index][i] +\
                                    data['ExEg']['matrix'][gamma_index][i+1])*0.5 \
                                    for (i,_) in enumerate(data['ExEg']['matrix']\
                                                        [gamma_index][:-1])]
data['ExEg']['exen centers'] = [(data['ExEg']['matrix'][exen_index][i] +\
                                    data['ExEg']['matrix'][exen_index][i+1])*0.5 \
                                    for (i,_) in enumerate(data['ExEg']['matrix']\
                                                        [exen_index][:-1])]
data['Singles'] = {'raw' : in_data['Singles'].to_numpy()} # axis 0 is the data, 1 is the bin edges
data['Singles']['histo'] =  data['Singles']['raw'][0] # axis 0 is the data, 1 is the bin edges
data['Singles']['centers'] = [(data['Singles']['raw'][1][i] +\
                                    data['Singles']['raw'][1][i+1])*0.5 \
                                    for (i,_) in enumerate(data['Singles']['raw'][1][:-1])]


# now we want to loop over the excitation energy matrix 
# and project it onto the gamma ray axis

"""
Assume:
    that the excitation energy is binned in sub-100 keV bins
    you want 100 keV bins
    you want binning centered at Sn

Make a list of two element lists that is like:
    [[lo_0, hi_0],...,[lo_N, hi_N]]
"""
step = 0.1
Sn = 9.15 # Mo96 Sn is 9.154 MeV, take midway point
# first put Sn at zero:
exen_basis = data['ExEg']['exeg axis'] - Sn

proj_idxs = []
c_val = min(exen_basis)
max_val = max(exen_basis)
centers = []
while c_val < max_val: # i don't like while loops, but here we go.
    n_val, n_idx = find_nearest(exen_basis, c_val)
    m_val, m_idx = find_nearest(exen_basis, c_val+step)
    c_val += step
    proj_idxs.append([n_idx, m_idx])
    centers.append((n_val + m_val)*0.5 + Sn)
proj_idxs = proj_idxs[:-1]
centers   = centers[:-1]
projs = {}
sings = {}

for i,idx_pair in enumerate(proj_idxs):
    # get these slices of the matrix:
    #m_slice = data['ExEg']['matrix'][data_index][idx_pair[0]:idx_pair[1]]
    # range is [)
    col_idxs = range(idx_pair[0], idx_pair[1]+1)
    if idx_pair[1]+1 == len(data['ExEg']['exeg axis']):
        col_idxs = range(idx_pair[0], idx_pair[1])
    
    m_slice = data['ExEg']['matrix'][data_index][:,col_idxs] # this is the gamma projection
    projs[centers[i]] = m_slice.sum(axis=1)

    # get the corresponding singles spectrum
    # add 1 because slices are [)
    sings[centers[i]] = data['Singles']['histo'][idx_pair[0]:idx_pair[1]+1]

fit_df = pd.read_csv('fit_gammas_95mo.csv')
fit_params = {}

# pars can vary between [(1.+lo)*p, (1+hi)*p], so if p = 0.05 we vary at the 5% level
ranges = {'A' : {'lo' : -0.6, 'hi' : 0.5},
          'sigma' : {'lo' : -0.5, 'hi' : 0.5},
          #'bg_slope' : {'lo' : -0.2, 'hi' : 0.2},
          #'bg_offset' : {'lo' : -0.2, 'hi' : 0.2},
          #'step_amp' : {'lo' : -0.2, 'hi' : 0.2},
          'mean' : {'lo' : -0.05, 'hi' : 0.05}, # should be well constrained, but don't be afraid to vary
          }

plot_save_dir  = 'saved_plots' # this is the directory in which the plots will be stored
if not os.path.isdir(plot_save_dir):
    os.mkdir(plot_save_dir)

plot_base_name = '96mo_pp_goddess' # this is the beginning of the filename for the plot
plot_save_flag = True # if this is true we save pngs of plots zoomed in on each fit.
fig,ax = plt.subplots(figsize=(10,6))

# we're going to read in an efficiency curve from a csv file
# the units should be the same as the units on your gamma spec, or
# you'll need to convert. The curve that's in the file now is a fake one that
# is all 10%. we're going to take the data and make an interpolated curve
# so we can evaluate the efficiency at any point.
# the csv file can have as many other columns in it as you like, but it needs to
# have the two columns called out below.
# this assumes fractional efficiency, but if you are wrong and like yours in %
# you just divide by 100.

eff_df = pd.read_csv('efficiency_curve.csv')
eff_func = scipy.interpolate.interp1d(eff_df['Egamma (keV)'], eff_df['efficiency'])

header = """# make_pgamma_plot.py - AR version Feb 23, 2023
          # Surrogate measurement of ('Mo96(p,p`))
          # proton beam energy 14.35 MeV
          # Gamma probability of GAM_EN_PH keV line in Mo96)
          # M. Williams analysis of surrogate measurement performed October 2022
          # using GODDESS located at ATLAS
          # Mo95 target isotopics (we'll fix this)
          # Mo92 = 0.28, Mo94 = 0.58, Mo95 = 96.8, Mo96 = 1.54, Mo97 = 0.36, Mo98 = 0.45, Mo100 < 0.1 all values percent
          # Mo96 Sn = 9.15434 MeV Sp = 9.2975 MeV Mo96(p,p) Q-value = 0 MeV
          # Zero excitation energy starts at a proton energy = XXXX
          # Zero neutron energy starts at an excitation energy = 9.15434 MeV
          #=====================================================
          # surGam Mo96 GAM_EN_PH keV Prob
          #=====================================================
          # This file was generated on the following date: DATE_PH
          #Energy dependant probability of this gamma ray
          #QUANTITY Probability
          #Eex,MeV  Prob,1  dProb,1 E_CN2,MeV
          #%.06f\\t%.06g\\t%.06g\\t%.06f
          """
out_dfs = {}
headers = {}
transition_info = {}
"""
tdf = fit_df.loc[fit_df['Egamma'] < 800 & fit_df['Egamma] > 600]
tdf = fit_df.loc[fit_df['background'].str.contains('linear)]

"""

for idx,row in fit_df.iterrows():
    for cen in projs.keys():
        # you really only need to fit near Sn, so I'm selecting a window so this 
        # is not too painful
        if cen > row['Ex Min'] and cen < row['Ex Max']: 
            #print('\n\n',row['Egamma'], row['fit lo'], row['fit hi'])
            # now we do the fits.
            # first get the spectrum in the fitting window:
            # we make it two bins bigger than we need...
            rng_lo, idx_lo = find_nearest(data['ExEg']['gamma centers'], row['fit lo'])
            rng_hi, idx_hi = find_nearest(data['ExEg']['gamma centers'], row['fit hi'])

            idxs = [idx_lo-2, idx_hi+2]
            while projs[cen][idxs[0]] == 0:
                # if there's nothing in the fit, wiggle it:
                idxs[0] = idxs[0]-1
                rng_lo = data['ExEg']['gamma centers'][idxs[0]]

            while projs[cen][idxs[1]] == 0:
                # if there's nothing in the fit, wiggle it:
                idxs[1] = idxs[1]+1
                rng_hi = data['ExEg']['gamma centers'][idxs[1]]


            # I think the default is [), and we want []
            xVals = np.asarray(data['ExEg']['gamma centers'][idxs[0]:idxs[1] + 1])
            yVals = np.asarray(projs[cen][idxs[0]:idxs[1] + 1])            
            #print(sum(yVals))
            yErrs = np.asarray(np.sqrt(yVals))
            # if there are zero entries in yErrs this gets unhappy, so we replace
            # zeros with the median error:
            m = np.median(yErrs[yErrs > 0])
            yErrs[yErrs == 0] = m

            fit_params[cen] = Parameters()
            est_params = estimate_parameters(xVals,yVals)            
            # add the estimated parameters
            for par in est_params.keys():

                fit_params[cen].add(par,value = est_params[par],
                             min = est_params[par]*(1. + ranges[par]['lo']),
                             max = est_params[par]*(1. + ranges[par]['hi']))
            
            #fit_params[cen].add('bg_slope', value = 0., vary = False)
            fit_params[cen].add('bg_offset', value = 10., vary = True)

            # add the mean
            fit_params[cen].add('mean',value = row['Egamma'],
                         min = row['Egamma']*(1. + ranges['mean']['lo']),
                         max = row['Egamma']*(1. + ranges['mean']['hi']))

            minner = Minimizer(fit_gamma,fit_params[cen],
                                fcn_args=(xVals,yVals,yErrs),
                                calc_covar=True)
            meth = 'leastsq'

            # calculate final result
            try:
                result = minner.minimize(method=meth)
                final = yVals + result.residual
                #print(result.params)
                print(fit_report(result))
            except Exception as X:
                # if there's a problem print the error and exit
                print(X)
                plt.cla()
                plt.step(data['ExEg']['gamma centers'], projs[cen],label='data')
                plt.step(xVals, gaussian(result.params,xVals),
                        label='gaus',
                        linestyle='dashed')
                #plt.step(xVals, linear_bg(result.params,xVals),
                #        label='linear bg',
                #        linestyle='dashed')

                #plt.step(xVals, step_bg(result.params,xVals),
                #        label='step bg',
                #        linestyle='dashed')


                ax.set_xlim([rng_lo, rng_hi])
                plt.show()
                exit()


            # now we clear the figure:
            plt.cla()
            ax.set_xlabel('E$_{\gamma}$ (keV)')
            ax.set_ylabel('Counts')
            
            #print(fit_report(result))
            plt.step(data['ExEg']['gamma centers'], projs[cen],label='data')
            plt.step(xVals, final, label='fit') # total fit
            plt.step(xVals, gaussian(result.params,xVals),
                        label='gaus',
                        linestyle='dashed')
            """
            #plt.step(xVals, linear_bg(result.params,xVals),
            #            label='linear bg',
            #            linestyle='dashed')
            
            plt.step(xVals, step_bg(result.params,xVals),
                        label='step bg',
                        linestyle='dashed')
            """
            ax.set_xlim([rng_lo-5,rng_hi+5])
            ax.legend(loc='best')
            filename = ('%s_exen_%.02f_egamma_%.02f'%(CN,cen,row['Egamma'])).replace('.','p') # get rid of decimal points
            plt.savefig(os.path.join(plot_save_dir,'%s.png'%(filename)),bbox_inches='tight')
            # now we integrate the gaussian and record the fit:
            total_coinc = sum(gaussian(result.params, xVals))
            # error can be statistical for now, but should include other uncertainty in the end
            # now we correct the volume for efficiency:
            eff_value = eff_func(row['Egamma'])
            d_total_coinc = np.sqrt(total_coinc)
            eff_cor_coinc = total_coinc/eff_value # do any other corrections (LT? you want to here)
            # now we can calculate the gamma emission probability
            # get the total singles in the bin:
            #total_singles = sum(data['Singles']['histo'][idxs[0]:idxs[1] + 1])
            total_singles = sum(sings[cen])
            d_total_singles = np.sqrt(total_singles)
            prob =  eff_cor_coinc/total_singles

            # so this is easier later... 
            dProb = prob*np.sqrt(pow(d_total_coinc/total_coinc,2.) + 
                                 pow(d_total_singles/total_singles,2.))
            print(total_coinc)
            print(d_total_coinc)
            print(prob)
            print(dProb)

            if row['Egamma'] not in out_dfs.keys():
                now = datetime.now()
                current_time = now.strftime('%Y-%m-%d %H:%M:%S')

                headers[row['Egamma']] = header.replace('GAM_EN_PH','%.02f'%(
                                row['Egamma'])).replace('DATE_PH',current_time)

                transition_info[row['Egamma']] = row['latex'] # this has the transition info
                out_dfs[row['Egamma']] = {'#Eex,MeV'  : [],
                                          'Prob,1'    : [],
                                          'dProb,1'   : [],
                                          'E_CN2,MeV' : []}

            out_dfs[row['Egamma']]['#Eex,MeV'].append(cen)
            out_dfs[row['Egamma']]['Prob,1'].append(prob)
            out_dfs[row['Egamma']]['dProb,1'].append(dProb)
            out_dfs[row['Egamma']]['E_CN2,MeV'].append(Sn - cen)


""" the following is for the output files Jutta uses for her fits. 
    if you use this format her life will be easier, which will be good 
    some of this stuff needs to be updated for 96Mo or 86Kr, but you can see
    how it works.
"""
base_file_name = '96Mo_pp_gam_prob_GAM_EN_PHkeV.dat'
prob_dir = 'fit_results'
if not os.path.isdir(prob_dir):
    os.mkdir(prob_dir)

for gam in out_dfs.keys():
    file_name = base_file_name.replace('GAM_EN_PH','%i'%(int(gam)))
    out = '%s\n'%(headers[gam])
    out_dfs[gam] = pd.DataFrame.from_dict(out_dfs[gam])
    print(out_dfs[gam].head())
    for idx,row in out_dfs[gam].iterrows():
        out = '%s%.06f\t%.06g\t%.06g\t%.06f\n'%(out,
                                                row['#Eex,MeV'],
                                                row['Prob,1'],
                                                row['dProb,1'],
                                                row['E_CN2,MeV'])


    with open(os.path.join(prob_dir,file_name),'w') as out_file:
        out_file.write(out)


plt.cla()
# now make the gamma emission prob plot:
for gam in out_dfs.keys():

    plt.errorbar(out_dfs[gam]['#Eex,MeV'], out_dfs[gam]['Prob,1'],
                    yerr=out_dfs[gam]['dProb,1'], xerr=step/2.,
                    label='$%s$'%(transition_info[gam]))

ax.legend(loc='best')
ax.set_yscale('log')

xlims = ax.get_xlim()
ylims = ax.get_ylim()

plt.plot([Sn]*2, ylims, linestyle='dashed')

ax.set_ylim(ylims)

ax.set_xlabel('Excitation Energy (MeV)')
ax.set_ylabel('Gamma-Ray Emission Probability')
plt.savefig('example.png',bbox_inches='tight')
plt.show()

