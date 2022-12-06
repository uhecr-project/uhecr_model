import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u
import h5py, os

from crpropa import *

N_grid = 50
Esrc_range = 10 ** np.linspace(1, 6, N_grid) * u.EeV
D_range = 10 ** np.linspace(0, np.log10(250), N_grid) * u.Mpc
Eth = 57 * u.EeV
N_sim = 1_000
sim_file_name = "output/crpropa_events_grid.txt"
particle_type = nucleusId(14, 7)
output_file_name = "output/crpropa_sim_table_N_maxEsrc_1e6.h5"

nitrogen_fraction = np.zeros((N_grid, N_grid))
source_particle_fraction = np.zeros((N_grid, N_grid))
mean_arrival_energy = np.zeros((N_grid, N_grid))
mean_arrival_energy_maxA = np.zeros((N_grid, N_grid))
mean_arrival_energy_mcA = np.zeros((N_grid, N_grid))
fraction_above_Eth = np.zeros((N_grid, N_grid))
mean_mass_number = np.zeros((N_grid, N_grid))
mean_energy_per_nucleon = np.zeros((N_grid, N_grid))


E_arr = []
for i, Esrc in enumerate(Esrc_range):
    for j, D in enumerate(D_range):

        sim = ModuleList()
        sim.add(SimplePropagation(1 * kpc, 50 * kpc))
        sim.add(Redshift())
        sim.add(PhotoPionProduction(CMB()))
        sim.add(PhotoPionProduction(IRB_Kneiske04()))
        sim.add(PhotoDisintegration(CMB()))
        sim.add(PhotoDisintegration(IRB_Kneiske04()))
        sim.add(NuclearDecay())
        sim.add(ElectronPairProduction(CMB()))
        sim.add(ElectronPairProduction(IRB_Kneiske04()))
        sim.add(MinimumEnergy(0.1 * EeV))

        obs = Observer()
        obs.add(ObserverPoint())
        output = TextOutput(sim_file_name, Output.Event1D)
        obs.onDetection(output)
        sim.add(obs)

        source = Source()
        source.add(SourcePosition(D.to_value(u.Mpc) * Mpc))
        source.add(SourceRedshift1D())
        source.add(SourceEnergy(Esrc.to_value(u.EeV) * EeV))
        source.add(SourceParticleType(particle_type))

        sim.setShowProgress(False)
        sim.run(source, N_sim, True)
        output.close()

        output.close()
        sim_data = np.genfromtxt(sim_file_name, names=True)

        Z = np.array([chargeNumber(int(id)) for id in sim_data["ID"].astype(int)])
        A = np.array([massNumber(int(id)) for id in sim_data["ID"].astype(int)])
        E = 10 ** (np.log10(sim_data["E"]) + 18)

        idx_max_A = A[A == max(A)]
        A_list = list(A)
        A_mc = max(set(A_list), key=A_list.count)
        idx_mc_A = A[A == A_mc]

        nitrogen_fraction[i][j] = len(A[(A == 14) & (Z == 7)]) / N_sim
        source_particle_fraction[i][j] = N_sim / len(A)
        mean_arrival_energy[i][j] = np.mean(E)
        mean_arrival_energy_maxA[i][j] = np.mean(E[idx_max_A])
        mean_arrival_energy_mcA[i][j] = np.mean(E[idx_mc_A])
        fraction_above_Eth[i][j] = len(E[E > Eth.to_value(u.eV)]) / len(E)

        mean_mass_number[i][j] = np.mean(A)
        mean_energy_per_nucleon[i][j] = np.mean(E / A)

        os.remove(sim_file_name)

with h5py.File(output_file_name, "w") as f:
    f.create_dataset("Esrc_range", data=Esrc_range.to_value(u.EeV))
    f.create_dataset("D_range", data=D_range.to_value(u.Mpc))
    f.create_dataset("Eth", data=Eth.to_value(u.EeV))
    f.create_dataset("nitrogen_fraction", data=nitrogen_fraction)
    f.create_dataset("source_particle_fraction", data=source_particle_fraction)
    f.create_dataset("mean_arrival_energy", data=mean_arrival_energy)
    f.create_dataset("mean_arrival_energy_maxA", data=mean_arrival_energy_maxA)
    f.create_dataset("mean_arrival_energy_mcA", data=mean_arrival_energy_mcA)
    f.create_dataset("fraction_above_Eth", data=fraction_above_Eth)
    f.create_dataset("mean_mass_number", data=mean_mass_number)
    f.create_dataset("mean_energy_per_nucleon", data=mean_energy_per_nucleon)
