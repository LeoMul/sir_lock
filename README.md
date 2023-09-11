# sir_lock
The implementation of lockdowns for the SIR model on small world networks; particularly for the large-deviation approach of obtaining the density of states.

This code is original but the idea serves as a continuation of the work published by Yannick Feld at https://journals.aps.org/pre/abstract/10.1103/PhysRevE.105.034313.
Yannick Feld should also be credited here, as he provided a lot of guidance in the development of the application.

A paper publishing the results of this work is in progress as of 09.09.23.

AS OF 9/09/23, the program is being updated such that the locked-down network is regenerated each time the system goes into lockdown. 
The following checklist will show the progress:

1. The simple-sampling functions ✓

2. The programs measuring the simple-sampling (should be one line changes) ✓ for SIMPLE SAMPLING, CRITICAL THRESH, LIFESPAN

3. Large deviations ✓ BY REMOVING THE LOCKDOWN MARKOV MOVES, and re-drawing the network at each lockdown.
