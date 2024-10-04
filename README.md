# iccee
Ice Climate Coupling for the Elevation of Exoplanets


*simu_ice.cmd* is the cycling job scheduler for beluga; it calls sequentially the GCM and a python script for ice accumulation.

There are 2 methods for the ice accumulation step:


## stepa (*change_elevation.py*)
Move forward to the next step when a fraction (e.g. 10%) of the total ice has accumulated.

Isostasy is computed in 2D.


## stept (*guess_elevation.py*)
Estimate the ice level based on the temperature vertical profile.

Isostasy is computed locally.
