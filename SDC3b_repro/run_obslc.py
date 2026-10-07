import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import py21cmfast as p21c
from py21cmfast import plotting
from astropy.cosmology import Planck15
from scipy import interpolate
import os, sys
import tools21cm as t2c
from astropy.io import fits
from tqdm import tqdm
import yaml
from py21cmfast import cache_tools

#venv
#/data-archive/sdc/venvs/t2c_21cmf

#----------------
# my functions
def mybin_lightcone_in_frequency(lightcone, z_low, box_size_mpc, dnu,output_frequencies):
    '''                                                                                                                  
    Bin a lightcone in frequency bins.                                                                                   
                                                                                                                         
    Parameters:                                                                                                          
        lightcone (numpy array): the lightcone in length units                                                           
        z_low (float): the lowest redshift of the lightcone                                                              
        box_size_mpc (float): the side of the lightcone in Mpc                                                           
        dnu (float): the width of the frequency bins in MHz                                                              
                                                                                                                         
    Returns:                                                                                                             
        * The lightcone, binned in frequencies with high frequencies first                                               
        * The frequencies along the line of sight in MHz                                                                 
    '''
    #Figure out dimensions and make output volume                                                                        
    cell_size = box_size_mpc/lightcone.shape[0]

    distances = t2c.z_to_cdist(z_low) + np.arange(lightcone.shape[2])*cell_size

    #print('distances before',distances)                                                                                 
    input_redshifts = t2c.cdist_to_z(distances)
    input_frequencies = t2c.z_to_nu(input_redshifts)
    #nu1 = input_frequencies[0]
    #nu2 = input_frequencies[-1]
    #output_frequencies = np.arange(nu1, nu2, -dnu)

    #sys.exit()                                                                                                          
    output_lightcone = np.zeros((lightcone.shape[0], lightcone.shape[1], \
                                 len(output_frequencies)))

    mapping=np.zeros(len(output_frequencies))
    
    for i in range(output_lightcone.shape[2]):
        nu = output_frequencies[i]
        idx = int(t2c.find_idx(input_frequencies, nu))
        output_lightcone[:,:,i] = lightcone[:,:,idx]
        mapping[i]=idx

    return output_lightcone,  mapping



def mybin_lightcone_in_mpc(lightcone, mapping,ncells):
# simply go back from the previous mapping
    
    output_lightcone=np.zeros((lightcone.shape[0], lightcone.shape[1],ncells))-100.

    print('mapping vs lightcone',mapping.shape,lightcone.shape[2])
    for i in range(lightcone.shape[2]):
        idx=int(mapping[i])
        try:
            output_lightcone[:,:,idx]=lightcone[:,:,i]
        except:
            print(idx,'out of bounds')
            
    # check for -100 flag value and fill with nearest slice
    for i in range(output_lightcone.shape[2]):
        if (output_lightcone[0,0,i]==-100):
            print('Warming: empty slice',i)
            output_lightcone[:,:,i]=output_lightcone[:,:,i-1]
    return output_lightcone


def myphysical_lightcone_to_observational(lightcone_freq, output_freqs,input_z_low, input_z_high, output_dnu, output_dtheta, fov_deg,input_box_size_mpc=None, verbose=True):

    if input_box_size_mpc == None:
        input_box_size_mpc = conv.LB

    #For each output redshift: average the corresponding slices                                                          
    t2c.print_msg('Making observational lightcone...')
    t2c.print_msg('Binning in frequency...')

    n_cells_theta = int(fov_deg*60./output_dtheta)                                                                  
    n_cells_nu = len(output_freqs)

    
    print('number of pixels in angular dimension',n_cells_theta)
    print('number of pixels in frequency dimension',n_cells_nu)


    t2c.print_msg('Binning in angle...')
    output_volume = np.zeros((n_cells_theta, n_cells_theta, n_cells_nu))

    
    for i in tqdm(range(n_cells_nu), disable=not verbose):
        if i%10 == 0:
            t2c.print_msg('Slice %d of %d' % (i, n_cells_nu))
        z = t2c.nu_to_z(output_freqs[i])
        
        output_volume[:,:,i] = t2c.physical_slice_to_angular(lightcone_freq[:,:,i], z, slice_size_mpc=input_box_size_mpc, fov_deg=fov_deg,dtheta=output_dtheta, order=2)

    return output_volume


def myobservational_lightcone_to_physical(observational_lightcone, input_freqs, input_dtheta, output_cell_size,verbose=True, order=2):
    '''                                                                                                                  
    Interpolate a lightcone volume measured in observational (angle/frequency)                                           
    units into  physical (length) units. The output resolution will be set                                               
    to the coarest one, as determined either by the angular or the frequency                                             
    resolution. The lightcone must have the LoS as the last index, with                                                  
    frequencies decreasing along the LoS.                                                                                
                                                                                                                         
    Parameters:                                                                                                          
        observational_lightcone (numpy array): the input lightcone volume                                                
        input_freqs (numpy array): the frequency in MHz of each slice along the                                          
            line of sight of the input                                                                                   
        input_dheta (float): the angular size of a cell in arcmin                                                        
        verbose (bool): show progress bar                                                                                
        order (int): The order of the spline interpolation, default is 2.                                                
            The order has to be in the range 0-5.                                                                        
            Use order=0 for ionization fraction data.                                                                    
                                                                                                                         
    Returns:                                                                                                             
        * The output volume                                                                                              
        * The redshifts along the LoS of the output                                                                      
        * The output cell size in Mpc                                                                                    
    '''
    assert input_freqs[0] > input_freqs[-1]
    assert observational_lightcone.shape[0] == observational_lightcone.shape[1]

    #Determine new cell size - set either by frequency or angle.                                                         
    #The FoV size in Mpc is set by the lowest redshift                                                                   

    dnu = input_freqs[0]-input_freqs[1]
    z_low = t2c.nu_to_z(input_freqs[0])
    fov_deg = observational_lightcone.shape[0]*input_dtheta/60.
    fov_mpc = fov_deg/t2c.angular_size_comoving(1., z_low) # this is the maximum box size in mpc that's common to all slices
    cell_size_perp = fov_mpc/observational_lightcone.shape[0]
    cell_size_par = t2c.nu_to_cdist(input_freqs[-1])-t2c.nu_to_cdist(input_freqs[-2])
    

    print('Making physical lightcone with cell size %.2f Mpc' % output_cell_size)
    #Go through each slice along frequency axis. Cut off excess and                                                      
    #interpolate down to correct resolution                                                                              
    n_cells_perp = int(fov_mpc/output_cell_size)
    output_volume_par = np.zeros((n_cells_perp, n_cells_perp, observational_lightcone.shape[2]))

    print(output_volume_par.shape)
    
    for i in tqdm(range(output_volume_par.shape[2]), disable=not verbose):
        z = t2c.nu_to_z(input_freqs[i])
        output_volume_par[:,:,i] = t2c.angular_slice_to_physical(observational_lightcone[:,:,i],\
                                                    z, slice_size_deg=fov_deg, output_cell_size=output_cell_size,\
                                                    output_size_mpc=fov_mpc, order=order)

    return output_volume_par


#-------------------------------------------------------
# MAIN ROUTINE STARTS
#________________________________________________________


paramfile = sys.argv[1]             # Name of the parameter file    

with open(paramfile,'r') as f:
    params = yaml.load(f, yaml.CSafeLoader)
    
#User parameters from the input file
fov_mpc = params['Grid']['boxsize']
npix_per = params['Grid']['meshsize']
print('fov mpc',fov_mpc)

input_dir= params['Output']['inputs_basename']
txt=input_dir.split("/")
run_name=txt[len(txt)-2]

#additional user parameters
# PRUDUCTION DS
fov_deg = 8. # in degrees
npix_deg=512 #sub-arcmin resolution  

#minicubes
#fov_deg=1.
#npix_deg=256



do_c2ray=False

#user settings
root_dir=input_dir#'/data-archive/sdc/SDC3/inference/'
path_output = root_dir#+'EoR_lightcones/'+run_name+'/'

files=['lc_dT_EoR_21cmf','lc_dT_EoR_C2R']

#frequency intervals 
fr_mins=[151,166,181] #min freq
fr_maxs=[166,181,196] #max_freq #original 15 MHz freq spacing

#fr_maxs=[161,176,191] #test PS for smaller freq range to limit evolution
#fr_maxs=[156,171,186] #test PS for smaller freq range to limit evolution


ra_deg = 0.     # in degrees
dec_deg = -30.  # in degrees

#--------------------------------------------

output_dtheta = (fov_deg /npix_deg) * 60.  
print('arcmin resolutions',output_dtheta)

min_freq_MHz=np.min(fr_mins)
max_freq_MHz=np.max(fr_maxs)


# k bin definition                                                                                              
kbins_fix_par=np.loadtxt('kpar_t2c_SDC3b.txt')
kbins_fix_per=np.loadtxt('kper_t2c_SDC3b.txt')

dkper=(kbins_fix_per[1]-kbins_fix_per[0])/2.
dkpar=(kbins_fix_par[1]-kbins_fix_par[0])/2.


kbins_fix_par=np.array(kbins_fix_par)-dkpar #t2c specifies the beginning of the bins, not centres              \
kbins_fix_per=np.array(kbins_fix_per)-dkper


#-------------------------------------------------------------------------------
print('setting cosmology')
# 21cmFAST default cosmology is Planck 2018 
# from https://arxiv.org/pdf/1807.06209.pdf
# Table 2, last column [TT+TE+EE+lowE+lensing+BAO]

    

# as in parameters.yml
Om0 = 0.30964
Ob0 = 0.04897              
H0 = 67.66


t2c.set_hubble_h(H0/100.)  
t2c.set_omega_matter(Om0)
t2c.set_omega_baryon(Ob0)


print('Adopted cosmology:')
print('Om',t2c.const.Omega0)
print('Ov',t2c.const.lam)
print('H0',t2c.const.H0)
print('h',t2c.const.h)
print(vars(p21c.CosmoParams))


for f in range(len(files)):

    tag=files[f]
    for ff in range(len(fr_mins)):

        min_freq_PS = fr_mins[ff]  # in MHz
        max_freq_PS = fr_maxs[ff] # in MHz
    

        cell_size_mpc=fov_mpc/npix_per

        #-------------------------------------------------------------------------------
        #compute cosmological quantities

        zmax=t2c.nu_to_z(min_freq_MHz)
        zmin=t2c.nu_to_z(max_freq_MHz)


        print('fov mpc',fov_mpc)

        fov_mpc_2 = t2c.deg_to_cdist(fov_deg, zmax)


        print('zmax=',zmax)
        print('fov_deg=',fov_deg)
        print('cdist=',t2c.z_to_cdist(zmax))
        print('fov_mpc check',fov_mpc,fov_mpc_2)
        
        print('lumdist',t2c.luminosity_distance(zmax))
        angular_size_deg = t2c.angular_size_comoving(fov_mpc, zmax)

        print('Fov that the cube fills',angular_size_deg)

        lc0 = np.load(path_output+tag+'.npy')


        print('lightcone shape',lc0.shape) 

        #here read the redshift list
        zs0 = np.loadtxt(path_output+'lc_redshifts.txt')

        print('cell_size_mpc',cell_size_mpc)
    
        #estimate suitable freq resolution to match
        dnu1=np.abs(t2c.z_to_nu(t2c.cdist_to_z(t2c.nu_to_cdist(min_freq_PS)+cell_size_mpc))-min_freq_PS)
        dnu2=np.abs(t2c.z_to_nu(t2c.cdist_to_z(t2c.nu_to_cdist(max_freq_PS)+cell_size_mpc))-max_freq_PS)

        print(' frequency resolution allowed by the moc resolution between ', dnu1, dnu2)
    
        dfreq_MHz=0.1 # this is fixed but it needs to be consistent with the estimated above

        # frequency range for the observed product I want
        freqs = np.arange(min_freq_PS, max_freq_PS-dfreq_MHz, dfreq_MHz)# the -dfreq_MHz is to avoid overlap between consecutive cubes and PS. the interval is intended including the start freq and excluding the end freq

    
        redshift = t2c.nu_to_z(freqs) #these are the redshifts corresponding to the equally spaced frequencies
        zmin_f = np.amin(redshift)
        zmax_f = np.amax(redshift)


        z_start_index=np.argmin(np.abs(zs0-zmin_f))
        z_end_index=np.argmin(np.abs(zs0-zmax_f))


        fmax_z=t2c.z_to_nu(zs0[z_start_index])
        fmin_z=t2c.z_to_nu(zs0[z_end_index])

        #here check if worth including a bit more
        print('the intended frequency range is',min_freq_PS,max_freq_PS-dfreq_MHz)
        print('the actual frequency range isolated using redshift is',fmin_z,fmax_z)

        print('Error in frequency',np.abs(fmin_z-min_freq_PS),np.abs(fmax_z-max_freq_PS))

        # select the relevant portion of the lighcone and redshift depending on the frequencies
        zs = zs0[z_start_index:z_end_index]
        lc = lc0[:,:,z_start_index:z_end_index]
        print('shape of lightcone before conversion',lc.shape) # 750 (512, 512, 750)
        print('min,max z',zs.min(), zs.max())

        # converting physical to observational coordinates - given cosmology is different here
        angular_size_deg = t2c.angular_size_comoving(fov_mpc, zs)
        print('Minimum angular size: {:.2f} degrees'.format(angular_size_deg.min()))
        print('Maximum angular size: {:.2f} degrees'.format(angular_size_deg.max()))
        
        physical_freq = t2c.z_to_nu(zs) # these are the frequencies corresponding to the native redshifts
        print('Minimum frequency gap in the physical light-cone data: {:.2f} MHz'.format(np.abs(np.gradient(physical_freq)).min()))
        print('Maximum frequency gap in the physical light-cone data: {:.2f} MHz'.format(np.abs(np.gradient(physical_freq)).max()))
        #Minimum frequency gap in the physical light-cone data: 0.12 MHz
        #Maximum frequency gap in the physical light-cone data: 0.20 MHz

        zmin = zs.min()
        zmax = zs.max()


        # compute the frequency vector from fmin,fmax,dfreq_MHz
        obs_freq = np.arange(max_freq_PS-dfreq_MHz, min_freq_PS-dfreq_MHz, -dfreq_MHz) # frequencies that I want. they need to be ordered in reverse to comply with z direction and exclude the last frequency which goes into the next cone
        
        # convert in the frequency dimension but remember the mapping
        lightcone_freq, mapping = mybin_lightcone_in_frequency(lc,zmin,fov_mpc,dfreq_MHz,obs_freq)

        # check the portion of the physical lightcone covered by the mapping, and cut if necessary
        i_start=int(np.min(mapping))
        i_end=int(np.max(mapping))

        print('cutting lightcone between',i_start,i_end,' of 0,',lc.shape[2]-1)
        lc=lc[:,:,i_start:i_end+1]        
        mapping=mapping-np.min(mapping) # if elements are eliminated at the start, mapping needs to reflect this




        
# here start commenting away if I want to exclude the derees transformation
#        '''
        print('run_phys2obs')
        # complete the transformation in the perp space
        obs_lc = myphysical_lightcone_to_observational(lightcone_freq,obs_freq,
                                                                 zmin,
                                                                 zmax,
                                                                 dfreq_MHz,
                                                                 output_dtheta,
                                                                 fov_deg,
                                                                 input_box_size_mpc=fov_mpc)

        print('done',obs_lc.shape)

        print('dimensions of the observational cube')
        print('number of frequency channels',len(obs_freq)) # 477
        print('min,max freq',obs_freq.min(), obs_freq.max()) # 80.69294013598301 199.692940135983
        print('shape',obs_lc.shape) # (512, 512, 477)
        print('FoV filled with real cube',t2c.angular_size_comoving(fov_mpc,zmax))
        #-------------------------------------------------------------------------------


        ## this is probably not required anymore
        ##the cube has an approximate frequency range wrt what required. Check and cut to the best size
        #print(np.min(obs_freq),np.max(obs_freq),min_freq_PS,max_freq_PS)
        ## obs_freq is the frequencies, compare to intended (freqs)
        #print('shape before',obs_lc.shape)
        #f_start_index=np.argmin(np.abs(obs_freq-max_freq_PS)) #frequency is reversed at this stage
        #f_end_index=f_start_index+freqs.shape[0]
        #print(f_start_index,f_end_index)
        #obs_lc=obs_lc[:,:,f_start_index:f_end_index]
        #obs_freq=obs_freq[f_start_index:f_end_index]
        #mapping=mapping[f_start_index:f_end_index]
        #print(np.min(obs_freq),np,max(obs_freq))
        #print('shape after',obs_lc.shape)
        #print(freqs.shape)
        #min_freq_PS_actual=np.min(obs_freq) #this difference is due to the redshifts
        #print('min_freq_PS_actual',min_freq_PS_actual)
        #
        # for output in fits format, rearrange axis direction and unit 

        # prepare the cone for saving in FITS
        lc_out = np.float32(obs_lc.transpose()[::-1])
        lc_out /= 1000.  # mK to K       
        
        xyrefpix = npix_deg/2.
        xypix_deg = fov_deg/npix_deg
        
        # cut to the desired number of pixels
        # if I process a subset of frequencies, a larger FoV than initially specified is retained
        # by cutting to HII_DIM I go back to the FoV I wanted
        
        npix_per_obs=obs_lc.shape[0]
        print('dimension of observed cube',npix_per_obs)
        print('wanted dimension',npix_deg)
        
        if (npix_per_obs>npix_per):
            #print('Warning: the processed cube was too big!!')
            lc_out=lc_out[:,int(npix_per_obs/2-npix_deg/2):int(npix_per_obs/2+npix_deg/2),int(npix_per_obs/2-npix_deg/2):int(npix_per_obs/2+npix_deg/2)]

        print('final shape for FITS cube',lc_out.shape)

        sdim = '%i' % int(npix_per)
        sfov = '%i' % int(fov_deg)
        
        outputname = path_output+'obs_'+tag+'_'+str(min_freq_PS)+"_"+str(max_freq_PS)+"_"+str(dfreq_MHz)+'.fits'
        
                
        hdu = fits.PrimaryHDU(lc_out)
        hdul = fits.HDUList([hdu])

        hdul[0].header.set('CTYPE1', 'RA---SIN')
        hdul[0].header.set('CTYPE2', 'DEC--SIN')
        hdul[0].header.set('CTYPE3', 'FREQ    ')
        hdul[0].header.set('CRVAL1', ra_deg)
        hdul[0].header.set('CRVAL2', dec_deg)
        hdul[0].header.set('CRVAL3', min_freq_PS*1.e6)# here use the approximate value or the actual frequency?
        hdul[0].header.set('CRPIX1', xyrefpix)
        hdul[0].header.set('CRPIX2', xyrefpix)
        hdul[0].header.set('CRPIX3', 1)
        hdul[0].header.set('CDELT1', -fov_deg/npix_deg)
        hdul[0].header.set('CDELT2', fov_deg/npix_deg)
        hdul[0].header.set('CDELT3', dfreq_MHz*1.e6)
        hdul[0].header.set('CUNIT1', 'deg     ')
        hdul[0].header.set('CUNIT2', 'deg     ')
        hdul[0].header.set('CUNIT3', 'Hz      ')
        hdul[0].header.set('BUNIT',  'K      ')
        
            
        hdul.writeto(outputname, overwrite=True)


        #---------------------------------------------------------------------
        # PS estimation on the observed cone
        #go back to comoving dimension                                                                                  


        print('transforming the lightcone into comoving units')
        # here I specify to go back to the original mpc cell size
        data=myobservational_lightcone_to_physical(obs_lc, obs_freq, output_dtheta,cell_size_mpc)
        
        # cutting the central part of the cube, to exclude the padding
        n_wanted=lc.shape[0] # since I am using the same cell size as before, I just need to select the same number of pixels
        npix_extra=data.shape[0]-n_wanted
        npix_extra_half=int(npix_extra/2.)
        #print(data.shape,n_wanted,npix_extra)
        #exit()
        if (npix_extra_half >0):
            print('cutting a padding of thickness',npix_extra)
            data=data[npix_extra_half:npix_extra_half+n_wanted,npix_extra_half:npix_extra_half+n_wanted,:]

        if (data.shape[0] != lc.shape[0]):
            print("Error: the final lightcone in mpc dimension is the worg shape")
            print(data.shape[0],lc.shape[0])
            #exit()

# '''
# here end commenting away of I want to exclude the degrees transformation
        #data=lightcone_freq

        print(np.min(mapping),np.max(mapping))
        
        n_cells=int(np.max(mapping)+1)

        data_phy2=mybin_lightcone_in_mpc(data,mapping,n_cells)
        cell_phy2=cell_size_mpc

        print('Comparison of mpc cube before and after transformation')
        print('dimension comparision',data_phy2.shape,lc.shape)
        print('Sums comparisoon',np.sum(data_phy2),np.sum(lc))


        # compute PS on cube on the same selection in freq
        
        #data_phy=lc[:,:,0:n_cells] #this should not be necessary anymore
        data_phy=lc
        cell_phy=cell_size_mpc
        min_freq_cdist = t2c.z_to_cdist(zmin)
        max_freq_cdist = t2c.z_to_cdist(zmax)
        # if I provide dimensions in Mpc/h the power spectra are in Mpc/h and so are the k bins
        box_perp=cell_phy*data_phy.shape[0]*t2c.const.h
        box_par=cell_phy*data_phy.shape[2]*t2c.const.h

        # power spectrum of the original cube in mpc
        Pk, kper_bins, kpar_bins = t2c.power_spectrum_2d(data_phy,kbins=[kbins_fix_per,kbins_fix_par],box_dims=[box_perp,box_perp,box_par],nu_axis=2,)#window='blackmanharris')          

        file = path_output+"Pk_"+tag+"_oncube"+str(min_freq_PS)+"_"+str(max_freq_PS)+"_"+str(dfreq_MHz)+".txt"
    
        np.savetxt(file, Pk*1.e-6, fmt='%e') #save in K


        box_perp=cell_phy2*data_phy2.shape[0]*t2c.const.h
        box_par=cell_phy2*data_phy2.shape[2]*t2c.const.h

        # power spectrum of the twice-transformed cube in mpc
        Pk2, kper_bins2, kpar_bins2 = t2c.power_spectrum_2d(data_phy2,kbins=[kbins_fix_per,kbins_fix_par],box_dims=[box_perp,box_perp,box_par],nu_axis=2,)#window='blackmanharris')                                                    
        print('done')
        file = path_output+"Pk_"+tag+"_oncone"+str(min_freq_PS)+"_"+str(max_freq_PS)+"_"+str(dfreq_MHz)+".txt"
        np.savetxt(file, Pk2*1.e-6, fmt='%e') #save in K
#        exit()

        #saving the mpc to frequency mapping information for optimal reproducibility
        file = path_output+"Mapping_"+tag+str(min_freq_PS)+"_"+str(max_freq_PS)+"_"+str(dfreq_MHz)+".txt"
        np.savetxt(file, mapping, fmt='%e')

    if (do_c2ray == False):
        exit()

