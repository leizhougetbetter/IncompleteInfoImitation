The code applies to Fig 2 and Fig 3 in the main text.

## Included files

#### Monte Carlo simulations
- *fpcCalculation.jl*: The main program of the Monte Carlo simulation used in the paper. One may use the code 'Julia fpcCalculation.jl' to run the program.

- *inputParameters.inp*: This input file contains the parameters for the simulation setting, including
  - networkType: the type of the network (RRG: random regular graph),
  - weightsDistType: the type of the edge weight distribution (homogeneous, uniform, or power-law),
  - numRepetitionInput: the repetition times for Monte Carlo simulations,
  - numGenerationInput: the total generation of a single simulation,
  - selectionIntensityInput: the selection intensity ($\delta$ in the main text),
  - numNeighChosenMain: the number of the amount of social information ($s$ in the main text), can be replaced by the file *sNumNeighobrInput.txt*,
  - SSEstimatedCritical: the critical value of benefit ($b$) for parwise social dilemmas or multiplication factor ($r$) for group social dilemma 
  - SSDistance: the distance from the start and end points to SSEstimatedCritical for simulations 
  - SSGap: the gap between points
  - TTvalue: the value of cost ($c$), 
  - theta: the relative importance of personal information ($\theta$), 

- *NetworkStructure_RRG_Weighted_Homogeneous.dlm*: the adjacency matrix of the network, the network is homogeneous and unweighted

- *regular_C=0.dlm*: the adjacency matrix of regular network with population size $N=10$ and clustering coefficient $C=0$.

- *regular_C=0.5.dlm*: the adjacency matrix of regular network with population size $N=10$ and clustering coefficient $C=0.5$.

- *sNumNeighobrInput.txt*: the number of the amount of social information $s$ of each node

#### Others
- *LICENSE*: MIT License

- *README.MD*: This file 


## Dependencies

Julia files was tested using Julia 1.8.4, Microsoft Visual Studio Code Version 1.77, and Microsoft Windows 11.

## Running the software

All files should be put in the same folder. Use 'Julia fpcCalculation.jl' to run the program

The program can perform simulations both under the pairwise and group social dilemmas. please set the kinds of social dilemmas on the line 170, 184, 270, and 280.

This program generates an output file named *ConditionalFixationProbTime_Z_Homogeneous_U* where *Z* represents the newtwork type and *U* the value of selection intensity $\delta$. 
Here is a demo for the output:
4.35216	0.00900374	450187	49549813	0
4.35216	0.01108856	554428	49445572	0
The items stands for: $b$ or $r$, fixation probability (The first row for cooperators $\rho_C$ and the second row for defectors $\rho_D$), invasionTimes, extinctionTimes, coexistenceTimes.

## License

See LICENSE for details.
