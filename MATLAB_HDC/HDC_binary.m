clear;
load('data/dataset.mat');

binary_functions;
%percentage of the dataset used to train the algorithm
learningFrac = 0.25;
%dimension of the hypervectors
D = 10000;
%number of classes to be classify
classes = 5; 
%percision used in quantization of input EMG signals
precision = 1;
%dimension of the N-grams
ngram = 10; 
%number of acquisition's channels
channels = 4;
%maximum level for the quantization 
maxL = 21; 

%Init Item Memories returns IM (to map all the elements in the system in 
%the HD space) and chAM (to map the signal level of the channels)  matrices.
%Set the variables "channels" and "maxL" respectively with the number of channels 
%and the maximum level of the signals (this level refers to the amplitude of 
%the envelope of the signals after the prepocessing, feature extraction and scaling).

[chAM, iMch] = initItemMemories (D, maxL, channels);

%downsample the dataset using the value contained in the variable "downSampRate"
downSampRate = 175;

[TS_COMPLETE_1, L_TS_COMPLETE_1] = downSampling (COMPLETE_1, LABEL_1, downSampRate);
[TS_COMPLETE_2, L_TS_COMPLETE_2] = downSampling (COMPLETE_2, LABEL_2, downSampRate);
[TS_COMPLETE_3, L_TS_COMPLETE_3] = downSampling (COMPLETE_3, LABEL_3, downSampRate);
[TS_COMPLETE_4, L_TS_COMPLETE_4] = downSampling (COMPLETE_4, LABEL_4, downSampRate);

%generate the training matrices using the learning rate contined in the
%variable "learningFrac"
[L_SAMPL_DATA_1, SAMPL_DATA_1] = genTrainData (TS_COMPLETE_1, L_TS_COMPLETE_1, learningFrac, '-------');
[L_SAMPL_DATA_2, SAMPL_DATA_2] = genTrainData (TS_COMPLETE_2, L_TS_COMPLETE_2, learningFrac, '-------');
[L_SAMPL_DATA_3, SAMPL_DATA_3] = genTrainData (TS_COMPLETE_3, L_TS_COMPLETE_3, learningFrac, '-------');
[L_SAMPL_DATA_4, SAMPL_DATA_4] = genTrainData (TS_COMPLETE_4, L_TS_COMPLETE_4, learningFrac, '-------');

downSampRate = 50;

[TS_COMPLETE_5, L_TS_COMPLETE_5] = downSampling (COMPLETE_5, LABEL_5, downSampRate);
[L_SAMPL_DATA_5, SAMPL_DATA_5] = genTrainData (TS_COMPLETE_5, L_TS_COMPLETE_5, learningFrac, '-------');
 
%This loop it is possible to test the classification capabilities with
%different N-grams. 
for N = 1 : ngram
    
    %Train the HDC algorithm using the matrices generated with the function
    %"genTrainData"
    [numpat_1, hdc_model_1] = hdctrain (L_SAMPL_DATA_1, SAMPL_DATA_1, chAM, iMch, D, N, precision, channels);
    [numpat_2, hdc_model_2] = hdctrain (L_SAMPL_DATA_2, SAMPL_DATA_2, chAM, iMch, D, N, precision, channels); 
    [numpat_3, hdc_model_3] = hdctrain (L_SAMPL_DATA_3, SAMPL_DATA_3, chAM, iMch, D, N, precision, channels);
    [numpat_4, hdc_model_4] = hdctrain (L_SAMPL_DATA_4, SAMPL_DATA_4, chAM, iMch, D, N, precision, channels); 
    [numpat_5, hdc_model_5] = hdctrain (L_SAMPL_DATA_5, SAMPL_DATA_5, chAM, iMch, D, N, precision, channels); 

   
    %Test the HDC model on the entire dataset. 
    [acc_ex, acc] = hdcpredict  (L_TS_COMPLETE_1, TS_COMPLETE_1, hdc_model_1, chAM, iMch, D, N, precision, classes, channels);
    accuracy(N,1) = acc;
 
    [acc_ex, acc] = hdcpredict  (L_TS_COMPLETE_2, TS_COMPLETE_2, hdc_model_2, chAM, iMch, D, N, precision, classes, channels);
    accuracy(N,2) = acc;
 
    [acc_ex, acc] = hdcpredict  (L_TS_COMPLETE_3, TS_COMPLETE_3, hdc_model_3, chAM, iMch, D, N, precision, classes, channels);
    accuracy(N,3) = acc;
 
    [acc_ex, acc] = hdcpredict  (L_TS_COMPLETE_4, TS_COMPLETE_4, hdc_model_4, chAM, iMch, D, N, precision, classes, channels);
    accuracy(N,4) = acc;
  
    [acc_ex, acc] = hdcpredict  (L_TS_COMPLETE_5, TS_COMPLETE_5, hdc_model_5, chAM, iMch, D, N, precision, classes, channels);
    accuracy(N,5) = acc;
   
end
 