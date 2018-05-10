# PULP-HD: Accelerating Brain-Inspired High-Dimensional Computing on a Parallel Ultra-Low Power Platform

This project provides an accelerator for HD Computing. The code is available for the execution on different architectures. High-level MATLAB code is used to validate the processing chain. The pre-processing functions are not provided and the input dataset is composed by the envelope extracted from the raw EMG data. Moreover, the C code is available for the execution on an ARM Cortex M4 (we tested the algorithm on an STM32F4-DISCOVERY). Another version of the application can be executed on the ultra-low power multi-core Wolf in sequential (using 1 core) and in multi-core (up to 8 cores), with and wothout built-ins to optimize the execution. The parallelization for the multi-core execution is done through OpenMP directives. 

## Results 
Some meaningful results are provided to show the capabilities of our implementation. 
PULPv3 with 4 cores achieves a 3.7x end-to-end speed-up and 2x energy saving compared to its single core execution. 
PULPv3 single core results slightly more complex than ARM Cortex M4 (1.2x), but, imposing a detection latency equal to 10ms, performance shows a power boost of 4.9x. This energy saving can be improved up to 9.9x using 4 cores, lowering the operational frequency and working at a voltage equal to 0.5V. Furthermore, the application was evaluated on Wolf obtaining a speed-up equal to 2.8x with 1 core and up to 18.4x using 8 cores. 

## Structure
All the functions that compose the application are commented with a brief descriptions and input/output arguments. 
### MATLAB
-**`MATLAB_HDC/data`**: contains the envelope derived from EMG signals acquired through 4 acquisition channels <br />
-**`MATLAB_HDC/binary_functions.m`**: contains all the functions needed to perform HD Computing. <br />
-**`MATLAB_HDC/HDC_binary.m`**: main script, which executes the processing starting from data acquired from 5 different subjects. <br />
-**`MATLAB_HDC/compress_hypervectors.m`**: this function compress a matrix or an array in 32-bit unsigned integer variables. <br />
-**`MATLAB_HDC/data_file_creator.m`**: with this script is possible to create the data.h file for the C implementation.   <br />
### STM32F4
-**`STM_HDC/inc/data.h`**: contains the IM and CIM matrices and a testing matrix.   <br />
-**`STM_HDC/inc/init.h`**: contains all the parameters. <br />
-**`STM_HDC/src/main.c`**: main function. <br />
-**`STM_HDC/src/aux_functions.c`**: here you can find all the function related to the HD computing algorithm.<br />
-**`STM_HDC/inc/aux_functions.h`**: definitions of the functions in “aux_functions.c”.<br />
-**`STM_HDC/src/associative_memory.c`**: function used to classify new samples.<br />
-**`STM_HDC/inc/associative_memory.h`**: definition of the function in “associative_memory.c”. <br />
### PULP
-**` PULP_HDC/inc/data.h`**: contains the IM and CIM matrices and a testing matrix.<br />
-**` PULP_HDC/inc/init.h`**: contains all the parameters.<br />
-**` PULP_HDC/src/main.c`**: main function.<br />
-**` PULP_HDC/src/aux_functions.c`**: here you can find all the function related to the HD computing algorithm.<br />
-**` PULP_HDC/inc/aux_functions.h`**: definitions of the functions in “aux_functions.c”.<br />
-**` PULP_HDC/src/associative_memory.c`**: function used to classify new samples.<br />
-**` PULP_HDC/inc/associative_memory.h`**: definition of the function in “associative_memory.c”.<br />
-**` PULP_HDC/Makefile`**: you can compile and run the application with "make clean all run" command. <br />

## Future Works
The next release will include an on-chip implementation of the training function and the pre-processing functions. In this way, it will be possible to train the algorithm and to classify new raw samples in real-time.   
STAY TUNED!

## Bugs Report and Clarifications. 
If you found some bugs or for any problems or questions, please contact us! 
For more information, you can read and (eventually) cite our paper: 

F. Montagna, A. Rahimi, S. Benatti, D. Rossi, L. Benini, “PULP-HD: Accelerating Brain-Inspired High-Dimensional Computing on a Parallel Ultra-Low Power Platform,” In IEEE/ACM Design Automation Conference (DAC), 2018. arXive preprint arXive: 1804.09123





