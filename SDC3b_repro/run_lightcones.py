# from the python notebook

#import libraries
import numpy as np, os, sys
import matplotlib.pyplot as plt
import py21cmfast as p21c
from tqdm import tqdm
import astropy.units as u
import astropy.constants as cst
import tools21cm as t2c
import yaml

paramfile = sys.argv[1]             # Name of the parameter file                                                                                        

with open(paramfile,'r') as f:
    params = yaml.load(f, yaml.CSafeLoader)

#User parameters from the input file                                                                                                                    
boxsize = params['Grid']['boxsize']
#npix_per = params['Grid']['meshsize']
input_dir= params['Output']['inputs_basename']
txt=input_dir.split("/")
run_name=txt[len(txt)-2]

print('boxsize=',boxsize)

# set inputs
#root_dir='/data-archive/sdc/SDC3/inference/'
root_dir=input_dir

#path_input = root_dir+'Full_outputs/'+run_name+'/'
path_input=root_dir

path_output = root_dir#+'EoR_lightcones/'+run_name+'/'

if not (os.path.exists(path_output+'/')):
    os.mkdir(path_output)


path_pyc2ray = path_input+'result_pyC2Ray/'
path_21cmfast = path_input+'result_21cmFAST/'

#todo: make it dependent on whether or not there is the C2ray folder
do_c2ray=False #c2ray is not always produced. 
z_min=6
z_max=10

print('C2ray',do_c2ray)

#i_plot = 32
#boxsize = 256 # cMpc #change this 

#redshift = np.loadtxt(path_input+'redshifts.txt')
redshift = np.loadtxt(path_input+'redshifts.txt')
#/data-hot/for-backup/sdc/SDC3/inference/3b_test1/outputs/redshifts.txt 

#21cmfast ionization fraction
redshift_21cmfast, glob_xHI_21cmfast = np.loadtxt('%sglob_xHI.txt' %path_input, unpack=True)
#c2ray ionization fraction


os.system('cp '+path_input+'redshifts.txt '+path_output)
os.system('cp '+path_input+'glob_xHI.txt '+path_output+'glob_xHI_21cmf.txt')


if (do_c2ray ==True):
    redshift_pyc2ray, tot_nHI, tot_phots, photion_pyc2ray, glob_xHII_pyc2ray = np.loadtxt('%sPhotonCounts2.txt' %path_pyc2ray, usecols=(0, 1, 2, 4, 6), unpack=True)
    glob_xHI_pyc2ray = 1 - glob_xHII_pyc2ray
    # save the C2ray neutral fraction in the same format as 21cmfast
    data = np.column_stack([redshift_pyc2ray,glob_xHI_pyc2ray])
    datafile_path = path_output+"glob_xHI_C2R.txt"
    np.savetxt(datafile_path , data, fmt=['%d','%d'])

    #print(redshift_pyc2ray)
    #print(redshift)
    #exit()
    
redshift = np.loadtxt(path_input+'redshifts.txt')
redshift=redshift[redshift >=z_min]
redshift=redshift[redshift <=z_max]
#print(redshift)
#exit()


glob_dT_21cmfast = np.zeros_like(redshift)
glob_dT_pyc2ray = np.zeros_like(redshift)
std_dT_21cmfast = np.zeros_like(redshift)
std_dT_pyc2ray = np.zeros_like(redshift)


print('`Making 21cmfast cone')


for i, z in enumerate(redshift):
    print('Analysing redshift',redshift[i])
    xHI_py21 = np.load('%sxHI_z%.3f.npy' %(path_21cmfast, z))
    
    ovrd_21cmfast = np.load('%sgrids/ovrdens_z%.3f.npy' %(path_input, z))
    
    dT_21cmfast = t2c.mean_dt(z)*xHI_py21*(1+ovrd_21cmfast) 

    np.save('%sdT_z%.3f.npy' %(path_21cmfast, z), dT_21cmfast)

    
    glob_dT_21cmfast[i] = dT_21cmfast.mean()

    
    std_dT_21cmfast[i] = dT_21cmfast.std()


files_21cmfast = np.array(['%sdT_z%.3f.npy' %(path_21cmfast, z) for z in redshift])
files_21cmfast_xh = np.array(['%sxHI_z%.3f.npy' %(path_21cmfast, z) for z in redshift])


print('Calling make lightcone')
lc_dT_21cmfast, lc_redshift = t2c.make_lightcone(files_21cmfast, 
                                                z_low=z_min, z_high=z_max, file_redshifts=redshift, 
                                                los_axis=2, interpolation='linear', 
                                                reading_function=np.load, box_length_mpc=boxsize)

print('xHI done')



lc_dT_21cmfast_xh, lc_redshift_dummy = t2c.make_lightcone(files_21cmfast_xh, 
                                                z_low=z_min, z_high=z_max, file_redshifts=redshift, 
                                                los_axis=2, interpolation='linear', 
                                                reading_function=np.load, 
                                                box_length_mpc=boxsize)


print('xH done')

# save and plot the lightcones
    
freqs = t2c.z_to_nu(lc_redshift)

np.savetxt(path_output+'lc_freqs.txt', freqs, fmt='%.3f', header='observed frequency [MHz]')
np.savetxt(path_output+'lc_redshifts.txt', lc_redshift, fmt='%.3f', header='lightcone redshift')

np.save(path_output+'lc_dT_EoR_21cmf.npy', lc_dT_21cmfast)
np.save(path_output+'lc_xH_EoR_21cmf.npy', lc_dT_21cmfast_xh)



if (do_c2ray==True):

    for i, z in enumerate(redshift):
        print('Analysing redshift',redshift[i])

        xHI_pyc = 1 - np.load('%sxfrac_%.3f.npy' %(path_pyc2ray, z))
    
        ovrd_21cmfast = np.load('%sgrids/ovrdens_z%.3f.npy' %(path_input, z))
        

        dT_pyc2ray = t2c.mean_dt(z)*xHI_pyc*(1+ovrd_21cmfast)

        np.save('%sdT_z%.3f.npy' %(path_pyc2ray, z), dT_pyc2ray)
    

        glob_dT_pyc2ray[i] = dT_pyc2ray.mean()
    

        std_dT_pyc2ray[i] = dT_pyc2ray.std()



    #Create lightcone with tools21cm

    files_pyc2ray = np.array(['%sdT_z%.3f.npy' %(path_pyc2ray, z) for z in redshift])

    files_pyc2ray_xh = np.array(['%sxfrac_%.3f.npy' %(path_pyc2ray, z) for z in redshift])



    lc_dT_pyc2ray, lc_redshift = t2c.make_lightcone(files_pyc2ray, 
                                                z_low=z_min, z_high=z_max, file_redshifts=redshift, 
                                                los_axis=2, interpolation='linear', 
                                                reading_function=np.load, 
                                                box_length_mpc=boxsize)





    lc_dT_pyc2ray_xh, lc_redshift_dummy = t2c.make_lightcone(files_pyc2ray_xh, 
                                                z_low=z_min, z_high=z_max, file_redshifts=redshift, 
                                                los_axis=2, interpolation='linear', 
                                                reading_function=np.load, 
                                                box_length_mpc=boxsize)






    np.save(path_output+'lc_dT_EoR_C2R.npy', lc_dT_pyc2ray)
    np.save(path_output+'lc_xH_EoR_C2R.npy', lc_dT_pyc2ray_xh)




