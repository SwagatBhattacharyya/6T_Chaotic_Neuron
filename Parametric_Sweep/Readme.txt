This folder contains two script files:
SweepNeuron_IF.m    :: Sweeps the input current to the neuron, computing the firing rate, Lyapunov exponents, and orbit diagram. This code corresponds to Fig. 7 in the paper.
SweepNeuron_Param.m :: Sweeps the MITE FG parameters: beta_m, beta_s, and V_FG0 computing the Lyapunov exponents and orbit diagram. You can change which variable is being swept by adjusting 'SweepType' in line 7. This code corresponds to Fig. 8 in the paper.

Please note that these code were tested on an Intel Core i9-13900KF (24-core) workstation with: 128 GB of DDR5 DRAM, Samsung 990 PRO SSD, and NVIDIA RTX 3060 Ti GPU. SweepNeuron_Param.m takes on the order of 10 minutes to execute on the aforementioned workstation, and SweepNeuron_IF takes on the order of 2 hours.

Please adjust the code as necessary to execute on your workstation in an acceptable timeframe. Note that the number of plotted points may need to adjusted for memory reasons if your workstation does not have a dedicated GPU. Execution time can be adjusted by changing the number of sweep points in:
SweepNeuron_IF:    InputCurrentLevels on line 45
SweepNeuron_Param: InputLevels on line 48, line 50, and line 52.