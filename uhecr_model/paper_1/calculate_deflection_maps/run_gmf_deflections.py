import argparse
import os
import pickle
import sys

import numpy as np
from scipy.optimize import root
from fancy import Data

import crpropa

known_ptypes = {
    "p": (1, 1),
    "He": (4, 2),
    "C": (12, 6),
    "N": (14, 7),
    "O": (16, 8),
    "Na": (24, 12),
    "Si": (28, 14),
    "Fe": (56, 26),
}


parser = argparse.ArgumentParser(
    description="Run fits to data a selection of catalog, "
    + "detector, mass, and analysis type."
)
parser.add_argument(
    "--detector",
    dest="detector_type",
    action="store",
    default="TA2015",
    type=str,
    help="The type of detector config (from TA2015, auger2014, auger2010)",
    choices=["TA2015", "auger2014", "auger2010"],
)
parser.add_argument(
    "--ptype",
    dest="ptype",
    action="store",
    default="p",
    type=str,
    help="Assumption for UHECR mass.",
    choices=list(known_ptypes.keys()),
)

parser.add_argument(
    "--gmf",
    dest="gmf",
    action="store",
    default="JF12",
    type=str,
    help="Selection of GMF model.",
    choices=["JF12", "TF17", "PT11"],
)

parser.add_argument(
    "--random_seed",
    dest="random_seed",
    action="store",
    default=19990308,
    type=int,
    help="Random seed for CRPropa.",
)


args = parser.parse_args()

"""Setting up"""
# Define file containing source catalogue information
uhecr_file = "../../data/UHECRdata.h5"

# make output directory if it doesnt exist
if not os.path.isdir("output"):
    os.mkdir("output")

detector_type = args.detector_type
assert detector_type in ["auger2010", "auger2014", "TA2015"]

# UHECR mass
ptype = args.ptype
assert ptype in known_ptypes

# GMF model
gmf = args.gmf
assert gmf in ["JF12", "PT11", "TF17"]

# # set random seed
random_seed = args.random_seed

print(detector_type, ptype, gmf, random_seed)


"""set detector and detector properties"""
if detector_type == "TA2015":
    from fancy.detector.TA2015 import Eth, M, alpha_T, detector_properties
elif detector_type == "auger2014":
    from fancy.detector.auger2014 import Eth, M, alpha_T, detector_properties
elif detector_type == "auger2010":
    from fancy.detector.auger2010 import Eth, M, alpha_T, detector_properties
else:
    raise Exception("Unknown detector type!")

# define output file
output_file = "output/" + f"defl_{detector_type}_{gmf}_{ptype}_{random_seed}.pkl"

# construct Dataset
data = Data()
data.add_uhecr(uhecr_file, detector_type, ptype=ptype)
data.add_detector(detector_properties)

# get uhecr detected arrival direction and detected energy
uhecr_coord = data.uhecr.coord
uhecr_energy = data.uhecr.energy

# get lon and lat arrays for future reference
# shift lons by 180. due to how its defined in mpl
uhecr_lons = np.pi - uhecr_coord.galactic.l.rad
uhecr_lats = uhecr_coord.galactic.b.rad

omega_true = np.zeros((len(uhecr_lons), 2))
omega_true[:, 0] = np.pi - uhecr_lons
omega_true[:, 1] = uhecr_lats

# number of uhecrs in the sample
N_uhecr = len(uhecr_coord)

# also get reconstruction uncertainties
arr_dir_unc = np.deg2rad(data.detector.coord_uncertainty)

# convert SkyCoord -> crpropa.Vector3d() to use with CRPropa
# coord is in galactic frame, lon and lat
# Angle definitions:
# CRPropa uses
#   longitude (phi) [-pi, pi] with 0 pointing in x-direction
#   colatitude (theta) [0, pi] with 0 pointing in z-direction
# matplotlib expects
#   longitude [-pi, pi] with 0 = 0 degrees
#   latitude [pi/2, -pi/2] with pi/2 = 90 degrees (north)
# SkyCoord uses:
#   lon: [0, 2pi]
#   lat: [-pi/2, pi/2]

uhecr_vector3d = []
for i, coord in enumerate(uhecr_coord):
    v = crpropa.Vector3d()
    v.setRThetaPhi(1, np.pi / 2.0 - coord.galactic.b.rad, np.pi - coord.galactic.l.rad)
    uhecr_vector3d.append(v)


# set up CRPropa simulation and initialize objects
sim = crpropa.ModuleList()

# setup magnetic field
if gmf == "JF12":
    gmf = crpropa.JF12Field()
    gmf.randomStriated(random_seed)
    gmf.randomTurbulent(random_seed)
elif gmf == "PT11":
    gmf = crpropa.PT11Field()
else:
    gmf = crpropa.TF17Field()

# Propagation model, parameters: (B-field model, target error, min step, max step)
sim.add(crpropa.PropagationCK(gmf, 1e-4, 0.1 * crpropa.parsec, 100 * crpropa.parsec))

# observer at galactic boundary (20 kpc)
obs = crpropa.Observer()
obs.add(crpropa.ObserverSurface(crpropa.Sphere(crpropa.Vector3d(0), 20 * crpropa.kpc)))
# obs.onDetection(TextOutput('galactic_backtracking.txt', Output.Event3D))
sim.add(obs)
print(sim)

# Set up composition
A, Z = known_ptypes[ptype]
pid = -crpropa.nucleusId(A, Z)

# CRPropa random number generator
crpropa_randgen = crpropa.Random(random_seed)

# Position of the Earth in galactic coordinates
pos_earth = crpropa.Vector3d(-8.5, 0, 0) * crpropa.kpc


# Some functions we need below
def get_time_delay(c):
    """Returns delay between entering the galactic disc
    and arrival at Earth through magnetic field."""
    return (
        (c.getTrajectoryLength() - c.current.getPosition().getDistanceTo(pos_earth))
        / crpropa.c_light
        / (60 * 60 * 24 * 365)
    )


def fischer_int(kappa, cos_thetaP):
    """Integral of vMF function over all angles"""
    return (1.0 - np.exp(-kappa * (1 - cos_thetaP))) / (1.0 - np.exp(-2.0 * kappa))


def fischer_int_eq_P(kappa, cos_thetaP, P):
    """Equation to find roots for"""
    return fischer_int(kappa, cos_thetaP) - P


# Number of samples to draw
Nsamples = 1000

# Random directions based on reconstruction uncertainty
omega_rand = np.zeros((N_uhecr, Nsamples, 2))
omega_gal = np.zeros((N_uhecr, Nsamples, 2))

# cos(theta), dot product
cos_thetas = np.zeros((N_uhecr, Nsamples))

# time delay in years between straight and deflected trajectory
time_delays = np.zeros((N_uhecr, Nsamples))

# Iterate over Nsamples and the detected events
for j in range(Nsamples):
    if gmf == "JF12" and j % 50 == 0:
        # Alternate seed for random jF12 component every 10 samples
        seed = np.random.randint(10000000)
        gmf.randomStriated(seed)
        gmf.randomTurbulent(seed)
    for i, arr_dir in enumerate(uhecr_vector3d):
        energy = (
            np.random.normal(
                loc=uhecr_energy[i],
                scale=data.detector.energy_uncertainty * uhecr_energy[i],
            )
            * crpropa.EeV
        )

        rand_arrdir = crpropa_randgen.randVectorAroundMean(arr_dir, arr_dir_unc)

        c = crpropa.Candidate(
            crpropa.ParticleState(pid, energy, pos_earth, rand_arrdir)
        )
        sim.run(c)

        defl_dir = c.current.getDirection()
        time_delays[i, j] = get_time_delay(c)

        # append longitudes and latitudes
        # need to append np.pi / 2 - theta for latitude
        # also append the randomized arrival direction in lons and lats
        omega_rand[i, j, :] = (
            rand_arrdir.getPhi(),
            np.pi / 2.0 - rand_arrdir.getTheta(),
        )
        omega_gal[i, j, :] = (
            defl_dir.getPhi(),
            np.pi / 2.0 - defl_dir.getTheta(),
        )
        # defl_lons[i, j] = defl_dir.getPhi()
        # defl_lats[i, j] = np.pi / 2.0 - defl_dir.getTheta()

        # evaluate dot product between arrival direction (randomized)
        # and deflected vector
        cos_theta = rand_arrdir.dot(defl_dir)
        cos_thetas[i, j] = cos_theta


# evaluate kappa_d from scalar product
# how this works is shown in solve_kappad.ipynb
kappa_gmf_rand = np.zeros((N_uhecr, Nsamples))

P = 0.683  # as defined in Soiaporn paper

for (i, j), cos_theta in np.ndenumerate(cos_thetas):
    sol = root(fischer_int_eq_P, x0=1, args=(cos_theta, P))
    # print(sol)   # check solution

    kappa_sol = sol.x[0]
    kappa_gmf_rand[i, j] = kappa_sol

    # print("kappa = ", kappa_sol)

# evaluate mean kappa for each uhecr
kappa_gmf = np.mean(kappa_gmf_rand, axis=1)


pickle.dump(
    (
        kappa_gmf,
        omega_gal,
        omega_rand,
        omega_true,
        kappa_gmf_rand,
        cos_thetas,
        time_delays,
    ),
    open(output_file, "wb"),
    protocol=-1,
)
