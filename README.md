# microswimmer-model
modelling distributions of microswimmers in 2D channel Poissieulle flows with chemotaxis.


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


