This project is about bifurcation analysis in a Wilson-Cowan model which was derived from finite-network biophysical model.

- network_simulation.m describes the biophysical finite-neuron network. compute_step_vec.m is called here to compute each individual step of the network.
- swipe_bt.m calls network_simulation.m for varying bt values to create plotsof peak frequency and GABA conductance over bt values.
- IN_neuron_ode_vec.m and PY_neuron_ode_vec.m contain the odes describing individual neurons.
- Compute_parameters.m is used to compute all wilson-cowan parameters apart from the weights. firing_rate.m and fit_sigmoid.m are helper functions that are called here.
- find_weights.m is used to find weight combinations that satisfy the derived inequalities.
  if 'mode' is set to 'GS' it performs grid search over different values of weights and creates a file with the resulting oscillating peak for each pair that results in oscillatory activity.
  if 'mode' is set ti 'IN' it runs the simulation once with the defined weights and plots the phase plane with trajectories for the specified initial point and the time plot for Inhibitory I(t) and excitatory E(t) activity.
- wilson_cowan.m contains the ode describing the Wilson - Cowan model.
- bifurcation_analysis.m performs bifurcation analysis on the tuned Wilson - Cowan model 
