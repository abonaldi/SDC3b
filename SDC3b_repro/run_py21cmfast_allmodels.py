#ANNA: pass the common params throgh yaml and the astro params directly. This is done to quicly create minicubes for all models


import numpy as np, os, sys
#import matplotlib.pyplot as plt
import py21cmfast as p21c
import h5py, pickle, yaml
from tqdm import tqdm
from scipy.stats import binned_statistic_dd

import astropy.units as u

# read pyc2ray parameter file
#paramfile = 'parameters.yml'
paramfile = sys.argv[1]             # Name of the parameter file    

row=sys.argv[2] #which row to read of C2ray_params

with open(paramfile,'r') as f:
    params = yaml.load(f, yaml.CSafeLoader)

#read the parameters for different runs from file and select a set
par_array=np.loadtxt('C2ray_params.txt',skiprows=1)
selected=par_array[int(row)-1,:]

#modify the params for this run


inc=int(selected[0]) #increment to random seed
tag=str(inc) #tag for the run

print('Now doing model',tag)


f0=selected[1]
g1=selected[2]
f0_esc=selected[3]
al_esc=selected[4]
alpha_h=selected[5]
mturn=selected[6]
min_halo=selected[7]
clumping=selected[8]
rmax=selected[9]

params['Sources']['f0']=float(f0)
params['Sources']['g1']=float(g1)
params['Sources']['f0_esc']=float(f0_esc)
params['Sources']['al_esc']=float(al_esc)
params['Sources']['alpha_h']=float(alpha_h)
params['Sources']['mturn']=float(mturn)
params['Sources']['minmass']=float(min_halo)
params['Sinks']['clumping']=float(clumping)
params['Sinks']['R_max_cMpc']=float(rmax)
params['Cosmology']['seed']=params['Cosmology']['seed']+inc
path_out = params['Output']['inputs_basename']+'model_'+tag+'/'
params['Output']['inputs_basename']=path_out
params['Output']['results_basename']=path_out+'result_pyC2Ray/'
params['Output']['sources_basename']=path_out+'sources/'
params['Output']['density_basename']=path_out+'grids/'

paramfile_new=paramfile+'_model'+tag
with open(paramfile_new, 'w') as outfile:
    yaml.dump(params, outfile, default_flow_style=False)

#exit()

# grid parameters
boxsize = params['Grid']['boxsize']
meshsize = params['Grid']['meshsize']
ntreads = params['Grid']['ntreads']

# star formation parameters
Nion = params['Sources']['Nion']
f_star10 = np.log10(params['Sources']['f0'])
#f_star10 = np.log10(f0)
alph_star = params['Sources']['g1']
#alph_star = g1

# escaping fraction
f0_esc = np.log10(params['Sources']['f0_esc'])
#f0_esc = np.log10(f0_esc)
alph_esc = params['Sources']['al_esc']
#alph_esc = al_esc

# this is the t_star that we fit to match the t_star of 21cmFAST
t_star = params['Sources']['alpha_h']
#t_star = alpha_h

# turnover mass (used only by 21cmfast)
mturn = params['Sources']['mturn']


# minimum halo mass  (used only by 21cmfast)
min_mhalo = float(params['Sources']['minmass'])
#min_mhalo = float(min_halo)

# constant cluping factor
clumping = params['Sinks']['clumping']
R_max_cMpc = params['Sinks']['R_max_cMpc']
#R_max_cMpc = rmax

# cosmological constants
Om0 = params['Cosmology']['Omega0']
Ob0 = params['Cosmology']['Omega_B']
h = params['Cosmology']['h']

#ANNA read seed from file
random_seed=params['Cosmology']['seed']

print('Grid:', params['Grid'])
print('Sources:', params['Sources'])
print('Sinks:', params['Sinks'])
print('Cosmology:', params['Cosmology'])

# Define output path
path_out = params['Output']['inputs_basename']

if not (os.path.exists(path_out)):
    print('---> Created directory: %s' %path_out)
    os.makedirs(path_out)
    os.makedirs(path_out+'grids/')
    os.makedirs(path_out+'sources/')
    os.makedirs(path_out+'result_21cmFAST/')
    os.makedirs(path_out+'img/') 
else:
    print('---> Directory already exist: %s' %path_out)


# user_paramser, cosmology and astrophysical parameters
cosmo_params = p21c.CosmoParams(OMb=Ob0, OMm=Om0, hlittle=h)
user_params = p21c.UserParams(HII_DIM=meshsize, BOX_LEN=boxsize, DIM=meshsize*3, USE_INTERPOLATION_TABLES=True, N_THREADS=ntreads, HMF=1, SAMPLER_MIN_MASS=min_mhalo, SAMPLER_BUFFER_FACTOR=1.2, AVG_BELOW_SAMPLER=False)
astro_params = p21c.AstroParams(HII_EFF_FACTOR=1, M_TURN=mturn, R_BUBBLE_MAX=R_max_cMpc, F_ESC10=f0_esc, ALPHA_ESC=alph_esc, F_STAR10=f_star10, ALPHA_STAR=alph_star, t_STAR=t_star, SIGMA_STAR=0.0)


#print(astro_params)
#exit()


# flag options
flag_options = p21c.FlagOptions(USE_HALO_FIELD=True, HALO_STOCHASTICITY=True, USE_MASS_DEPENDENT_ZETA=True, USE_MINI_HALOS=False, INHOMO_RECO=True, 
                                USE_UPPER_STELLAR_TURNOVER=False, USE_EXP_FILTER=True, CELL_RECOMB=True)





# Save parameters to a file
with open('%scosmo_par.pkl' %path_out, 'wb') as file:
    pickle.dump(cosmo_params.defining_dict, file)

# Save parameters to a file
with open('%sastro_par.pkl' %path_out, 'wb') as file:
    pickle.dump(astro_params.defining_dict, file)

with open('%suser_par.pkl' %path_out, 'wb') as file:
    pickle.dump(user_params.defining_dict, file)

with open('%sflag_opt.pkl' %path_out, 'wb') as file:
    pickle.dump(flag_options.defining_dict, file)

# change the cache path
#random_seed = 918
path_cache = path_out+'21cmFAST-cache/'
p21c.config['direc'] = path_cache

lc_filename=path_out+'result_21cmFAST/lightcone.save'


# Results are better for the halo finder using the top-hat filter for reionisation
p21c.global_params.HII_FILTER = 0

# set clumping factor
p21c.global_params.CLUMPING_FACTOR = clumping

# set ionizing photons per baryons
p21c.global_params.Pop2_ion = Nion

print('running initial conditions')
# run initial conditions
initial_conditions = p21c.initial_conditions(user_params=user_params, cosmo_params=cosmo_params, random_seed=random_seed)
print('done')

print('running lightcone')
# run 21cmfast lightcone
lc = p21c.run_lightcone(redshift=5., init_box=initial_conditions, flag_options=flag_options, astro_params=astro_params, lightcone_quantities=('xH_box', 'density', 'brightness_temp'), regenerate=False, random_seed=random_seed, direc=path_cache)
print('done')

print('saving to file')
#ANNA: here save the lightcone to compute the PS later
lc.save(lc_filename)
print('lightcone saved')


# plot lightcone
#plt.figure(figsize=(18, 5))
#plt.pcolormesh(lc.lightcone_redshifts, np.linspace(0, 256, 128), lc.xH_box[32], cmap='jet')
#plt.xlim(5, 14.4)
#plt.savefig(path_out+'img/slice_21cmfast.png', bbox_inches='tight')

# read and save cache files to be used with pyc2ray
redshifts = []
nr_halo_arr = []
print('saving info for C2Ray')
for i in tqdm(range(lc.global_xH.size)):
    print('index',i)
    # get redshift of the chace
    z = lc.cache_files['ionized_box'][i][0]

    with h5py.File(lc.cache_files['halobox'][i][1]) as f:
        n_ion = f['HaloBox']['n_ion'][:]
    #np.save('%ssources/nion_z%.3f.npy' %(path_out, z), n_ion)

    # number of halos
    nr_halo = np.count_nonzero(n_ion)
    nr_halo_arr.append(nr_halo)
    
    # save number of ionizing photons for pyC2Ray
    if(nr_halo != 0):
        # redshift of the sources
        redshifts.append(z)

        # store only indexes and Nion
        ii, jj, kk = np.nonzero(n_ion)
        Nion = n_ion[ii, jj, kk]

        np.savetxt('%ssources/src_z%.3f.txt' %(path_out, z), np.array([ii, jj, kk, Nion]).T, fmt='%d\t%d\t%d\t%.2f', header='i\tj\tk\tn_ion = Mstar/Vcell*Nion*fesc [Msun/Mpc^3]')
    else:
        pass

    # access the chace file
    with h5py.File(lc.cache_files['pt_halos'][i][1]) as f:
        halo_mass = f['PerturbHaloField']['halo_masses'][:]
        halo_pos = f['PerturbHaloField']['halo_coords'][:]
        nr_halo = halo_mass.size
        
    # save number of ionizing photons for pyC2Ray
    if(nr_halo > 0): 
        np.savetxt('%ssources/halo_z%.3f.txt' %(path_out, z), np.hstack((halo_pos, halo_mass[..., None])), fmt='%d\t%d\t%d\t%.2f', header='i\tj\tk\thalo mass [Msun]')
    else:
        pass

    # save Overdensity field
    with h5py.File(lc.cache_files['perturb_field'][i][1]) as f:
        overdens = f['PerturbedField']['density'][:].astype('float64')

    np.save('%sgrids/ovrdens_z%.3f.npy' %(path_out, z), overdens)

    # save Ionization field
    with h5py.File(lc.cache_files['ionized_box'][i][1]) as f:
        xHI = f['IonizedBox']['xH_box'][:]
    np.save('%sresult_21cmFAST/xHI_z%.3f.npy' %(path_out, z), xHI)
    
print('done')
print('saving final info')
# store redshift 
all_redshifts = np.array([lc.cache_files['ionized_box'][i][0] for i in range(len(lc.cache_files['halobox']))])
redshifts = np.array(redshifts) #all_redshifts.copy()
nr_halo_arr = np.array(nr_halo_arr)
np.savetxt('%sall_redshifts.txt' %path_out, all_redshifts, fmt='%.3f')
np.savetxt('%sredshifts.txt' %path_out, redshifts, fmt='%.3f')
np.savetxt('%sglob_xHI.txt' %path_out, np.array([all_redshifts, lc.global_xH]).T, fmt='%.3f\t%.5e')
print('done')
