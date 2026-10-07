import numpy as np, sys, os
import time
import pyc2ray as pc2r
import yaml

from c2ray_sdc3b import C2Ray_SDC3b
# ======================================================================
# Test run for pyc2ray for the SDC3b
# ======================================================================

# Global parameters
num_steps_between_slices = 2        # Number of timesteps between redshift slices
paramfile = sys.argv[1]             # Name of the parameter file

row=sys.argv[2] #which row to read of C2ray_params

par_array=np.loadtxt('C2ray_params.txt',skiprows=1)
selected=par_array[int(row)-1,:]

inc=int(selected[0]) 
tag=str(inc) #tag for the run
paramfile=paramfile+'_model'+tag

#modify the redshift zred_0 keyword
with open(paramfile,'r') as f:
    params = yaml.load(f, yaml.CSafeLoader)

path=params['Output']['inputs_basename']
redshifts_21cmfast=np.loadtxt(path+'redshifts.txt')
print(redshifts_21cmfast[0])
params['Cosmology']['zred_0']=float(redshifts_21cmfast[0])

#update the paramfile with the new zred_0
with open(paramfile, 'w') as outfile:
    yaml.dump(params, outfile, default_flow_style=False)

#print(paramfile)
#exit()
# this runs the simulation
sim = C2Ray_SDC3b(paramfile=paramfile)

#modify the redshift
#print(sim.zred_0)
#exit()

# copy parameter file in output directory
os.system('cp %s %s' %(paramfile, sim.results_basename))

# Get redshift list (test case)
zred_array = np.loadtxt(sim.inputs_basename+'redshifts.txt', dtype=float)

# check for resume simulation
if(sim.resume):
    i_start = np.argmin(np.abs(zred_array - sim.zred))
    sim.resume = i_start+1
else:
    i_start = 0

# Measure time
tinit = time.time()

# Loop over redshifts
for k in range(i_start, len(zred_array)-1):
    zi = zred_array[k]       # Start redshift
    zf = zred_array[k+1]     # End redshift

    pc2r.printlog(f"\n=================================", sim.logfile)
    pc2r.printlog(f"Doing redshift {zi:.3f} to {zf:.3f}", sim.logfile)
    pc2r.printlog(f"=================================\n", sim.logfile)

    # Compute timestep of current redshift slice
    dt = sim.set_timestep(zi, zf, num_steps_between_slices)

    # Read input files
    sim.read_density(fbase='ovrdens_z%.3f.npy' %zi, z=zi)

    # Read source files
    #srcpos, normflux = sim.ionizing_flux(fbase='halo_z%.3f.txt', z=zf, dt=dt)
    #srcpos, normflux = sim.read_sources_diff(fbase='nion_z%.3f.npy', zi=zi, zf=zf, dt=dt*num_steps_between_slices)
    #srcpos, normflux = sim.read_sources(fbase='nion_z%.3f.npy' %zf, z=zf)
    srcpos, normflux = sim.read_sources(fbase='src_z%.3f.txt' %zf, z=zf)
    
    # Write output
    sim.write_output(z=zi, ext='.npy')

    # Set redshift to current slice redshift
    sim.zred = zi

    # Loop over timesteps
    for t in range(num_steps_between_slices):
        tnow = time.time()
        pc2r.printlog(f"\n --- Timestep {t+1:n}. Redshift: z = {sim.zred : .3f} Wall clock time: {tnow - tinit : .3f} seconds --- \n", sim.logfile)

        # Evolve Cosmology: increment redshift and scale physical quantities (density, proper cell size, etc.)
        sim.cosmo_evolve(dt)

        # Evolve the simulation: raytrace -> photoionization rates -> chemistry -> until convergence
        sim.evolve3D(dt, normflux, srcpos)

    # Evolve cosmology over final half time step to reach the correct time for next slice (see note in c2ray_base.py)
    sim.cosmo_evolve_to_now()

# Write final output
sim.write_output(zf, ext='.npy')
tend = time.time()

pc2r.printlog(f"\n --- End of the simulation. Wall clock time: {tend - tinit: .3f} seconds --- \n", sim.logfile)
