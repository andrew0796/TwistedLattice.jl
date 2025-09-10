## SU(3) Cooling
Run code with `julia --threads=n SU3_Cooling.jl` where `n` is the number of threads you want to use

Output is saved into an HDF5 file

## General
This is a much more complete working example, run with `julia --threads=n general_cooling.jl` where `n` is the number of threads you want to use. This will create two data files which are cooled from the same initial configuration using two cooling methods (SU(2) subgroups and polar decomposition). Also created is text files which have the progress bars printed out to them. To view the progress in real time, in another terminal run `tail -f SU3_24_6_6_6_twists_14_23_{METHOD}_progress.txt` where `{METHOD}` should be replaced by either `SU2_subgroups` or `polar_decomposition`. Note that this is run in order, where `SU2_subgroups` is run first.

## Higher Charge Instantons
One can use the improved action to isolate instantons with topological charge greater than the minimal charge as described in [hep-lat/9309009](https://arxiv.org/abs/hep-lat/9309009). Run this using `julia --threads=n higher_charge_cooling.jl` where `n` is the number of threads you want to use

## Progress Bars
Both of these examples will produce progress bars which show the number of iterations, the action decomposed into electric + magnetic parts, and the log of the absolute value of the convergence of the action.