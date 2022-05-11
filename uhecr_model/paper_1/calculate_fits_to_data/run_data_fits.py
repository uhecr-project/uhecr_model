import os
import argparse
from fancy import Data, Model, Analysis


parser = argparse.ArgumentParser(
    description="Run fits to data a selection of catalog, "
    + "detector, mass, and analysis type."
)
parser.add_argument(
    "--source",
    dest="source_type",
    action="store",
    default="SBG_23",
    type=str,
    help="The source catalogue used (from SBG_23, 2FHL_250Mpc, swift_BAT_213)",
    choices=["SBG_23", "2FHL_250Mpc", "swift_BAT_213"],
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
    "--model",
    dest="model_type",
    action="store",
    default="joint",
    type=str,
    help="The stan model considered for simulation (from joint, joint_gmf)",
    choices=["arrival", "joint", "joint_gmf"],
)
parser.add_argument(
    "--ptype",
    dest="ptype",
    action="store",
    default="p",
    type=str,
    help="Assumption for UHECR mass.",
    choices=["p", "N", "Fe"],
)

parser.add_argument(
    "--gmf",
    dest="gmf",
    action="store",
    default="JF12",
    type=str,
    help="Selection of GMF model.",
    choices=["JF12", "TF17", "PT11", None],
)

parser.add_argument(
    "--random_seed",
    dest="random_seed",
    action="store",
    default=19990308,
    type=int,
    help="Random seed for Stan.",
    # choices=["JF12", "TF17", "PT11"],
)


args = parser.parse_args()

# parser.add_argument(
#     "--check",
#     dest="check",
#     action="store_true",
#     help="Check if table is correctly constructed.",
# )


"""Setting up"""
# Define location of Stan files
stan_path = "../../stan/"

# Define file containing source catalogue information
source_file = "../../data/sourcedata.h5"
uhecr_file = "../../data/UHECRdata.h5"

# make output directory if it doesnt exist
if not os.path.isdir("output"):
    os.mkdir("output")

# source_types = ["SBG_23", "2FHL_250Mpc", "swift_BAT_213"]
source_type = args.source_type

# detector_types = ["auger2010", "auger2014", "TA2015"]
# detector_type = "auger2014"
detector_type = args.detector_type

# UHECR mass
ptype = args.ptype

# model (arrival, joint, joint + GMF)
model_type = args.model_type

# GMF model
gmf = args.gmf

# set random seed
random_seed = args.random_seed

# flag to control showing plots or not
show_plot = False

print(source_type, detector_type, ptype, model_type, gmf)


"""set detector and detector properties"""
if detector_type == "TA2015":
    from fancy.detector.TA2015 import detector_properties, alpha_T, M, Eth
elif detector_type == "auger2014":
    from fancy.detector.auger2014 import detector_properties, alpha_T, M, Eth
elif detector_type == "auger2010":
    from fancy.detector.auger2010 import detector_properties, alpha_T, M, Eth
else:
    raise Exception("Unknown detector type!")

# define output file
output_file = (
    "output/"
    + f"{model_type}_fit_{source_type}_{detector_type}_{gmf}_{ptype}_{random_seed}.h5"
)

# construct Dataset
data = Data()
data.add_source(source_file, source_type)
data.add_uhecr(uhecr_file, detector_type, ptype=ptype)
data.add_detector(detector_properties)
# data.show();

# table file
table_file = "../../tables/tables_{0}_{1}.h5".format(source_type, detector_type)

if model_type == "arrival":
    """Fit arrival model"""
    model_fname = stan_path + "arrival_direction_model.stan"
    summary = b"Fit of the arrival direction model to data"
elif model_type == "joint":
    model_fname = stan_path + "joint_model_tightB.stan"
    summary = b"Fit of the joint model to data"
elif model_type == "joint_gmf":
    model_fname = stan_path + "joint_gmf_model_tightB.stan"
    summary = b"Fit of the joint + GMF model to data"
else:
    raise Exception(f"Undefined model type {model_type}")

# construct arrival model obejct
model = Model(model_filename=model_fname, include_paths=stan_path)
model.compile()
model.input(Eth=Eth)  # EeV

# Define an Analysis object to bring together Data and Model objects
analysis = Analysis(
    data,
    model,
    analysis_type=model_type,
    filename=output_file,
    summary=summary,
)

# Define location of pre-computed values used in fits
# (see relevant notebook for how to make these files)
# Each catalogue has a file of pre-computed values
analysis.use_tables(table_file)

# Fit the Stan model
fit = analysis.fit_model(chains=8, iterations=500, seed=random_seed)

# Save to analysis file
analysis.save()
