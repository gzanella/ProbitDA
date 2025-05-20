This folder contains the code to replicate Tables 1 and 2 in Analysis of data-augmentation samplers for high-dimensional probit regression (Ascolani and Zanella).

It includes the following files:

- "Functions_cpupling.R", which contains the functions to implement all the couplings.

- "Simulations_gprior.R", which contains the code to replicate the first two rows of Tables 1 and 2.

- "Simulations_noint.R", which contains the code to replicate the third and fourth row of Tables 1 and 2.

- "Simulations_withint.R", which contains the code to replicate the last row of Tables 1 and 2.


See Section 6 for the details of the simulation.

NOTE: all the files admit the possibility of generating imbalanced data (Table 1) or data generated from the model (Table 2). In the code this is highlighted by the comment "#generate the y: COMMENT THE ONE YOU DO NOT USE"