# microswimmer-models
modelling distributions of microswimmers in 2D channel 
Poissieulle flows with chemotaxis.

### IBM_noRT_base_model
this is the base model that is mostly taken from Smitha's
paper. This is used to build the other models by adding the
tumbling. Here the parameters are used as the base parameters
as well in order to have as a reference point.

### memory_model




Normal model with no run and tumble:
- swimmer_trajectory_DPP
can now store multiple trajectories for 
a y by theta grid of initial points
also produces the period histogram which
may make the period_histogram obsolete

- Channel_crossing
uses NaN, or averages from values that are not NaN
finds first hit time/angle of each wall
finds crossing channel time

wall hitting angle distributions

- period_histogram  (for reorientation times)


- first_hit_time (stop iteration after first hit)



- diffusion measure


- angle_capped_model

- angle_uncapped_model
i believe this does not work as this chooses a delay time,
then it does a normal run and then flips at the delay time,
then this 

i have tried different parameters
and attempted it to force

- angle_instantaneous

- memory_alt_params

- memory_model

- continuous_memory
this is the main model where i suppose that the levels of
the phosphorolated proteins dissipate exponentially

