# Notes on how this subset has been compiled from the repo

Findings:

- Figure 3, 11, 12 in uhecr_model/figures/verification/figures_kappa_gmf.ipynb
- Figure 1 in uhecr_model/figures/verification/figures_TA_SBG_data.ipynb
- Figure 2, powerpoint

# Folders
## `calculate_fits_to_data`

> Fits all models to data using STAN. All tables must be prepared for this.

A more compact copy of notebooks/gmf/fit_to_data with scripts to mass-produce the fits for all combinations of models and assumptions on a condor cluster.

ToDo:

- [X] Fix parallelization problem -> Fixed at local condor level
- [X] Are 8 chains enough?
- [ ] There is a figures folder that need to go somewhere else once Figures are separated from computations
- [X] Do not run nuclei for arrival only model (because just scales B/Z)
- [ ] Generate and use CRPropa loss tables
- [ ] It's not clear which GMF model is selected for stan and how it happens. We should do JF12 and TF17.

## `calculate_deflection_maps`

> Generates deflection maps and compute event catalogs of kappa_d values. 

Collection of scripts for the backpropagation in CRPropa3 based on the previous [eval_kappad notebook](../notebooks/gmf/eval_kappad.ipynb).

Several things need to be improved:

- [X] There is no sampling of energy uncertainty to generate the maps, needs to be added according to detector resolution
- [X] Directions are sampled from a gaussian instead of Fisher, replace with randFisherVector from CRPropa -> No need since error it is quoted exp. error and assumed to be Gaussian
- [X] Need to include TF17 field
- [X] Calculate time delays (curiosity, just a note to myself)
- [X] Increase number of samples and parallelize
- [ ] Omega true needs to be subtracted from np.pi for plotting. Maybe fix it here but then needs to be done consistently everywhere.

## `deflection_figures_exposure

Plots that have been used in the intro of the paper and for explaining some details.

- [ ] Not sorted or cleaned up, yet :)

## `calculate_loss_function`

Location for new scripts to calculate an interpolated energy loss function using CRPropa3.

Items:

- [ ] interpolate $E_{\rm Earth}(E_\text{source}, D_\text{source}, A_\text{source})$
- [ ] cross check $\langle A_{\rm Earth}\rangle(E_\text{source}, D_\text{source}, A_\text{source})$ or $\langle \ln{A} \rangle$

## `simulations_for_appendix`

Here, we need some visualization of the simulation of the GMF case that used galactic lenses.