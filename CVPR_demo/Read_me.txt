%% 
Run Ward_demo1.m to generate results on WARD dataset. Generated results will be
slightly different(better!) compared to the results reported in CVPR. 

%%Description of variables as outputs in workspace. 
Data_set= Input dataset. 
Train_single structure : single shot data. cam cell contains features of persons from three cameras. camp contains person identitites. 
Train_multi structure: Training data in multi-shot format. cam contains 3  by 2 cells, corresponding to three camera pairs (1-2,1-3,2-3) ,i.e.,
each row corresponds to each pair of cameras. camp contains identities from three camera pairs. 
Test_multi sturucture : Test data. same format as Train_multi. 

 



MA_store= 4x1 cell for 4 methods (Exact, Greedy, Half and Baseline), each cell contains
three binary matrices(in sparse format). Three matrices corresponds to Manual labeling for three camera pairs (1-2,1-3,2-3) and (i,j)-th element corresponds to the pair consists of
i-th person in camera 1and j-th person in camera 2.  Matrx element value 1 denotes manual label needs to be collected for that specific pair; 0 otherwise. 

FL_store= 4x1 cell for 4 methods (Exact, Greedy, Half and Baseline), each cell contains
three matrices(in sparse format). Three matrices corresponds to Aggregated Labeling Matrices (that is these matrices which labels we have access to
after label propagation) for three camera pairs (1-2,1-3,2-3). Matrix element value 1/0 denotes label of that specfic pair and  NaN indicates we do not have
access to the label of that specific pair.  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
Dependencies: 
1) CPLEX software must be installed. 
2) KISSME toolbox[1] (already included). 


If you use our code or algorithm, please cite the following paper:

Roy S, Paul S, Young NE, Roy-Chowdhury AK. "Exploiting Transitivity for Learning Person Re-identification Models on a Budget." 
In Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition 2018 (pp. 7064-7072).


contact: sroy004@ucr.edu

References:
[1] M. Koestinger, M. Hirzer, P. Wohlhart, P. M. Roth, and H. Bischof. "Large scale metric learning from equivalence constraints". In CVPR, 2012