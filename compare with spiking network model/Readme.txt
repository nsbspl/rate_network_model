This folder contains the code simulating the spiking network model (STAR Methods S6), whose results were compared with our rate network model (Figure 7). Users need to install NEST to run the code in Python. 

Users can simulate the results for different Vim-DBS frequencies by using the corresponding DBS-induced post-synaptic current (I_DBS) into Vim neurons. I_DBS in response to the DBS frequencies {5, 10, 20, 30, 50, 100, 200Hz} are provided in this folder. Users can simulate I_DBS from arbitrary DBS frequencies with the codes in the folder "supporting codes/Tsodyks & Markram model of STP".

In line 18, users can provide the specific DBS frequency to simulate. In line 131, users need to adjust the appropriate directory for his/her own computer.