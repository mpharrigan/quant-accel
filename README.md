MAccelerator: Quantitative Analysis of Adaptive Sampling for MSMs
==================================================================

MSMs model the dynamics of a system by discritizing conformational space
into a finite number of states and counting the number of transitions between
each state. These raw counts are used to estimate reversible transition
probabilities. Theoretically, these raw counts can also be used to estimate
uncertainty in the model and select optimal states from which we can spawn
additional simulation. In practice, many design choices can affect the
speed-up gained from adaptive sampling.

## Tunable parameters for investigation
 - :white_check_mark: Adaptive frequency -- length of each trajectory before
                      re-starting simulation
 - :white_check_mark: Degree of Parallelization -- Number of trajectories
                      to run simultaneously
 - :white_circle:     Adaptive lag-time -- The system has a natural lag-time
                      for converged models. Would it help to use a different
                      lag-time for generating the counts matrix from which
                      we adapt?

## Adaptive schemes
 - :white_check_mark: Uniform -- Select a new state at random.
 - :white_circle:     Min-Counts -- Sort states by counts and choose states
                      in order
 - :white_circle:     Weighted-Counts -- Select new states with a probability
                      as counts^(-n) where n parameterizes exploration vs.
                      refinement

## Toy systems
 - :white_check_mark: Simulating from a known transition matrix
    - 20 state alanine dipeptide
 - :white_circle:     OpenMM simulation
    - Muller potential


# Software Design

This package takes a very object-oriented approach in its design. Each
combination of [toy system, adaptive scheme] is enclosed in a 
`Configuration` object. Each `Configuration` has a

    - Method that yields a number of `AdaptiveParam`
      objects which contain tunable parameters.
    - `Simulator` object
    - `Modeller` object
    - `ConvergenceChecker` object
    - `Adapter` object

## Installing

Running
```
python setup.py install
```
will install all relevant code, as well as the code needed to generate
the reference data for the included toy systems. Use
```
make_reference_data.sh
```
to build the toy systems. You may need to run `setup.py` again to
copy the reference data to the install location.

## Running
the script `maccel.py` can be used via command line to generate

    - A sample python configuration script
    - A sample job script to be used for `qsub`

The python script controlls any final details of a system before running.
For example, it makes sense to define the "grid" of tunable parameters
in this configuraiton script. The adaptive scheme will probably be
specified here as well, although right now it is just taken to be the default.
