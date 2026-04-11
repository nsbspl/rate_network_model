Steps for running the route optimization to replicate the results in paper:

1. The inputs are: 
(1) The experimental (reference) firing rate computed by time histogram: FR_time_hist_conc.mat
(2) The DBS-induced post-synaptic current, obtained by the Tsodyks & Markram model: I_DBS_conc.mat
(3) The initial parameters obtained from anatomical literature and our previous model: ini_para_network_model.mat
Users can see some more details of explanation of the inputs in the Readme file in the "input data" folder.

2. Run Stabilizing_Function_preliminary_fit.mat with the above 3 inputs to get the model parameters 
from the preliminary fit.

3. Use the parameters from Step 2 as input parameters to run route_opt_Global_Stage.mat.

4. Use the final parameters from Step 3 as input parameters to run route_opt_Refining_Stage.mat.

Additional note: 

a). To replicate the results presented in Table S4, the users do 4 iterations in Step 3, and then 
move to Step 4 for the remaining iterations. 

b). The results from the route optimization are robust to the variabilities from simulating our previous model and
the slight changes (e.g., number of iterations, proportion of high frequency DBS data) during the route optimization process.
