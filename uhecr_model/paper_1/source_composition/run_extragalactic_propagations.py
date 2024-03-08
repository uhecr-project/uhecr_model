import argparse
import os
import pickle

import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u
import h5py

from crpropa import *

known_ptypes = {
    "p": (1, 1),
    "He": (4, 2),
    "C": (12, 6),
    "N": (14, 7),
    "O": (16, 8),
    "Ne": (20, 10),
    "Na": (24, 12),
    "Al": (26, 13),
    "Mg": (26, 13),
    "Si": (28, 14),
    "Cl": (35, 17),
    "Ca": (40, 20),
    "Cr": (52, 24),
    "Fe": (56, 26),
}

parser = argparse.ArgumentParser(
    description="Run forward simulations of UHECR propagations" +
    "for a given source composition on a rigidity-distance grid."
)

parser.add_argument(
    "--ptype",
    dest="ptype",
    action="store",
    default="Si",
    type=str,
    help="Assumption for UHECR mass.",
    choices=list(known_ptypes.keys()),
)

args = parser.parse_args()

'''Setting up'''
# grid size
N_grid = 30

# Range of source rigidities
Rsrc_range = 10**np.linspace(1, 4, N_grid) * u.EV

# Range of source distances
D_range = 10**np.linspace(0, np.log10(100), N_grid) * u.Mpc

# number of particles to simulate
N_sim = 1_000

# doesnt matter since it will access it for each iteration
sim_file_name = f"output/crpropa_events_grid_{args.ptype}.txt"

# arrays to store data into
mean_charge = np.zeros((N_grid, N_grid, len(known_ptypes)))
mean_mass_number = np.zeros((N_grid, N_grid, len(known_ptypes)))
mean_energy_per_nucleon = np.zeros((N_grid, N_grid, len(known_ptypes)))
mean_rigidity = np.zeros((N_grid, N_grid, len(known_ptypes)))
mean_rigidity_nuclei = np.zeros((N_grid, N_grid, len(known_ptypes)))
mean_A_Aarr = np.zeros((N_grid, N_grid, len(known_ptypes)))

# define globally the arrival particle A and Z
A_arr, Z_arr = known_ptypes[args.ptype]

# define array of all particle types that we need to forward propagate
source_particles = {}

for (pt, ptup) in list(known_ptypes.items()):
    # all particles such that mass number is larger than the mass number of 
    # the arrival composition
    if ptup[0] > A_arr:
        source_particles[pt] = ptup

# construct equivalencies to convert between rigidity and energy
rigidity_energy =[
    (u.V, u.eV, lambda x : x, lambda x : x)
]

# run the simulation for each source composition

for k, (pt, ptup) in enumerate(list(source_particles.items())):
    for i, Rsrc in enumerate(Rsrc_range):
        print(f"Particle: {pt}, Rsrc: {Rsrc:.2f}")
        for j, D in enumerate(D_range):
            
            # construct simulation
            sim = ModuleList()
            sim.add( SimplePropagation(1*kpc, min([D.to_value(u.Mpc)/10, 10]) * Mpc) )
            sim.add( Redshift() )
            sim.add( PhotoPionProduction(CMB()) )
            sim.add( PhotoPionProduction(IRB_Kneiske04()) )
            sim.add( PhotoDisintegration(CMB()) )
            sim.add( PhotoDisintegration(IRB_Kneiske04()) )
            sim.add( NuclearDecay() )
            sim.add( ElectronPairProduction(CMB()) )
            sim.add( ElectronPairProduction(IRB_Kneiske04()) )
            sim.add( MinimumEnergy( 1 * EeV) )

            # add observer and output
            obs = Observer()
            obs.add( ObserverPoint() )
            output = TextOutput(sim_file_name, Output.Event1D)
            obs.onDetection( output )
            sim.add( obs )

            # define source properties (source comp., position, energy)
            source = Source()
            source.add( SourcePosition(D.to_value(u.Mpc) * Mpc) )
            source.add( SourceRedshift1D() )
            particle_type = nucleusId(*ptup) 

            # should set Z * R (grid should be in R), then use energy here
            Esrc = Rsrc * ptup[1]  # E = R * Z

            source.add( SourceEnergy(Esrc.to_value(u.EeV, equivalencies=rigidity_energy) * EeV) )
            source.add( SourceParticleType(particle_type) ) 
            
            # run simulation
            sim.setShowProgress(False)
            sim.run(source, N_sim, True)
            output.close()

            # now open simulation file
            sim_data = np.genfromtxt(sim_file_name, names=True)

            # get observational parameters
            Z = np.array([chargeNumber(int(id)) for id in sim_data['ID'].astype(int)])  
            A = np.array([massNumber(int(id)) for id in sim_data['ID'].astype(int)])  
            E = 10**(np.log10(sim_data['E']) + 18)

            # want R such that all Z != 0
            Z_nonzero_indces = np.argwhere(Z != 0)
            R = E[Z_nonzero_indces] / Z[Z_nonzero_indces]  # rigidity
            E_per_A = E / A

            # append all physical variables

            # average arrival rigidity, charge, mass number, energy per nucleon
            mean_charge[i][j][k] = np.mean(Z)
            mean_mass_number[i][j][k] = np.mean(A)
            mean_energy_per_nucleon[i][j][k] = np.mean(E/A)
            mean_rigidity[i][j][k] = np.mean(R)
            
            # mean rigidity for those that are not protons A > 1
            mean_rigidity_nuclei[i][j][k] = np.mean(R[A[Z_nonzero_indces] > 1])

            # number of particles that are same as arrival composition that we want
            mean_A_Aarr[i, j, k] = len(Z == Z_arr)


# store data in hdf5 file

# get mass and charge number separately
source_particles_A = list(ptup[0] for ptup in list(source_particles.values()))
source_particles_Z = list(ptup[1] for ptup in list(source_particles.values()))
# Save
with h5py.File(f"output/extragal_propa_{args.ptype}.h5", "w") as f:
    f.create_dataset("src_comp_labels", data=list(source_particles.keys()))
    f.create_dataset("src_comp_A", data=source_particles_A)
    f.create_dataset("src_comp_Z", data=source_particles_Z)
    f.create_dataset("Rsrc_range", data=Rsrc_range.to_value(u.EV))
    f.create_dataset("D_range", data=D_range.to_value(u.Mpc))
    f.create_dataset("mean_mass_number", data=mean_mass_number)
    f.create_dataset("mean_charge", data=mean_charge)
    f.create_dataset("mean_energy_per_nucleon", data=mean_energy_per_nucleon)
    f.create_dataset("mean_rigidity", data=mean_rigidity)
    f.create_dataset("mean_rigidity_nuclei", data=mean_rigidity_nuclei)
    f.create_dataset("mean_A_Aarr", data=mean_A_Aarr)