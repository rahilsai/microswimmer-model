# Modelling Microswimmers in 2D channel poissieulle flows
modelling distributions of microswimmers in 2D channel 
Poissieulle flows with chemotaxis. In this model I am taking a linear
chemotactic gradient. I start with a base model of comparison, then I 
consider different methods of introducing the run and tumbles. Finally 
I use diffent methods and plots to compare the distributions and measure 
effects of the run and tumble.

## Models
### IBM_noRT_base_model (need to plot)
this is the base model that is mostly taken from Smitha's
paper. This is used to build the other models by adding the
tumbling. Here the parameters are used as the base parameters
as well in order to have as a reference point.

### angle_instantaneous (single plot of y-theta)
The tumble rate is calculated based on the angle, so when facing
the upper wall the sin value is the greatest resulting in a 
minimum rate of tumbling and vice versa for facing the lower wall.
This causes a bias towards the upper wall.

Then for the actual tumbling, I take the probability of flipping 
at each step as the rate multiplied by the dt as when the dt is 
small enough I then get the average number of steps within 
P(flip|theta) = lambda(XX(2)) * dt

### memory_model_delay (capped)
Here I use the T as the length of a period then divide it by 
nSteps to get the dt (timesteps). At the start of each period
I calculate a delay jump time tau based on the concentrations 
measured at at the current and previous period. I then divide 
this time by dt to calculate the number of steps needed and use 
a ceilingfunction so that the flips occur at a specific timestep. 
If the number of delay step are greater than the number of steps 
in a period (so if the tau>T), the jump just never happens and 
a new concentration or attractant is measured.

### angle_capped_model
Similarly to the memory,I use T as the length of a period and so
this is the maximum delay time. Here an angle is sampled at the 
start of each period and then the tumble rate is calculated from
this and is then used to calculate the delay time.

### angle_uncapped_model
The same as angle capped model except yhere is no limit to the 
value the delay time tau can be.

I believe this does not work as this chooses a delay time,
then it does a normal run and then flips at the delay time,
then this 

I have tried different parameters
and attempted it to force

- continuous_memory
this is the main model where i suppose that the levels of
the phosphorolated proteins dissipate exponentially

## Plots/figure methods
Normal model with no run and tumble:

- diffusion measure

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

