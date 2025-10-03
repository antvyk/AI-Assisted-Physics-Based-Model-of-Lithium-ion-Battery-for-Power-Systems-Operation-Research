# AI-Assisted-Physics-Based-Model-of-Lithium-ion-Battery-for-Power-Systems-Operation-Research
This repository contains code and supporting materials for the paper AI-Assisted Physics-Based Model of Lithium-ion Battery for Power Systems Operation Research. The code implements a mixed-integer linear programming reformulation of neural-network model of lithium-ion battery energy storage that includes battery operation and degradation. The neural-network model is trained on datasets generated from a physics-based digital twin. 

The energy arbitrage application of LiBESS is used to demonstrate the capabilities of the proposed model in optimization. This GitHub repository contains: (1) code to derive the strategic operation of LiBESS using the proposed model, (2) code to simulate lithium-ion cell behavior with a single-particle model coupled with solid electrolyte interphase growth, (3) code to train neural networks, and (4) code to generate the training dataset. To run the optimization problem, the following parameter files are used:
- Alberta’s electricity market pool prices in 2022
- The system level parameters of LIBESS
- The cell-level parameters

## Requirements:
The code for both simulation and optimization models is written in the Julia programming environment. To run an optimization problem, a valid Gurobi/CPLEX license is required.
