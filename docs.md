# Guide to using this package

## Step-by-step instructions

1. Create the environment using the .yml file.
2. Install both `stan_utility` and `fancy` from the `uhecr-project` organization (via pip)


**Note:** skip steps 3, 4 if the necessary table and UHECR datasets are already present (and no new settings are created (ex. if you want to perform the analysis with a new dataset, you need to run these steps)).


3. Run `uhecr_model/tables/precompute_tables.py` to precompute the tables for kappa and energy evaluations.
4. Run `data/precompute_kappa_gmf.ipynb` to precompute kappa_GMF for each UHECR dataset.
5. Initialize the runs by `run.py`. Here is how to use this:

    a. Modify the configurations in line 104 for custom runs
    Customizable variables are:
    - sources
    - detectors
    - sim_models (models used to simulate)
    - sim_models_for_fit (model of simulation that fit_models will use)
    - fit_models (models that will fit the simulation / data)
    - ptypes (particle type)
    - seeds
    - end_labels (additional config to modify filename)

    b. Run in command line with possible flags
    - mode: controls the type of run you want to do
        - sim: simulation
        - fit_sim: fit the simulation
        - fit_data: fit the data
        - fit: fit both simulation and data
        - sim_and_fit: simulate and fit the simulation
        - all: do it all
    - run_all: perform run with max. configuration
    - debug: debug mode
    - yes: skip confirmation
    - dryrun
    
        This should launch runs that will distribute to different cores (based on what its doing). You can check the progress in `outputs/progress_report.txt`

6. Once the runs are complete, plot the results from `plot_results.py`. 
    a. Configure at line 113 in the same fashion as above. Only difference is that we can add modify the header as a label for the folders.

    b. Again, run this from the command line. Flags are:
    - mode: the type of plots to construct
        - sim, fit_sim, fit_data as above
        - defl: for deflection skymaps
        - all
    - dist_mode: type of comparison to make with f-distribution
        - model, source, detector, ptype, seed
    - plot_all: plot all plots from all possibe configurations
    - debug
    - dryrun
    - yes

    This should launch some commandline text indicating which plot is constructed until it finishes.

    This would create a `plots` directory in the parent directory that will divide the plots based on the header file you gave it. The directory structure basically goes like this:
```
header
    - data
        - arrival_direction
        - joint
        - joint_gmf
        - fdist
    - simulation
        - joint
            - arrival_direction
            - joint
            - joint_gmf
            - fdist
        - joint_gmf 
            - arrival_direction
            - joint
            - joint_gmf
            - fdist
```

Note: the second tree after simulation is the **simulation**, and the tree after that is the **fits**.


## What does the code do under this?

### For `runs.py`

- uses `fit_model.py` and `simulate_data.py`, which is based on the underlying notebooks from the original uhecr_model from Fran.

#### For **simulation**

1. Sets the initial parameters (detector, alpha_T, M, Eth, L, F0) based on the source / detector
2. Compiles stan file containing the model we want to simulate
3. simulate the data from stan.
4. saves to simulation output file

#### For **fitting**

1. Initialize parameters from simulation output file (for simulations) or data file (for data)
2. Compile the model from the stan file
3. set up the analysis and use the precomputed tables
4. fit the model with the simulation / data
5. save the file into output file 

## Some ending remarks
- the GMF sims are not integrated with this workflow. Instead, use the notebook `figures/simulate_and_plot_defl.ipynb` and edit the configurations over there to construct them.

## Some future improvements
- some naming for the directories can be improved as they can be misleading. This should be changed to be more concise.
- the progress checking could be done in a better way (ex. what if we want to run one configuration, then run another one while it is still running?)
- some notebooks are still messy. These can and should be cleaned up.
- Some good documentation should be constructed (probably as a Jupyter notebook?)




Let me know if there are any other details I should provide. Thanks!

