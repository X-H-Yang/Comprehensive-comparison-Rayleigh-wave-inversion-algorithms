function [Y_hat,all_index,NumFea_hat,NumWin_hat,NumEnhan_hat] = ...
    bls_regression_Y_sub_noTest(train_x,train_y,Fea_vec,Win_vec,...
    Enhan_vec,validation_x,Vp_true,den_true,f,modes_num_vec,index_vec,...
    index_vec_all,mean_y,std_y,validation_x_norm)
%
% note: built upon previously released codes (https://broadlearning.ai/)
%
% Chen and Liu 2017: Chen, C. P., & Liu, Z. (2017). Broad learning system: 
%                    An effective and efficient incremental learning system
%                    without the need for deep architecture. IEEE 
%                    transactions on neural networks and learning systems, 
%                    29(1), 10-24.
%
% Chen et al. 2018: Chen, C. P., Liu, Z., & Feng, S. (2018). Universal 
%                   approximation capability of broad learning system and 
%                   its structural variations. IEEE transactions on neural 
%                   networks and learning systems, 30(4), 1191-1204.
%
assert(isfloat(train_x), 'train_x must be a float');

%% training
C = 2^-30;      %----C: the regularization parameter for sparse regualarization
s = .8;              %----s: the shrinkage parameter for enhancement nodes
best = 1;
result = [];

MAXGEN = length(Fea_vec)*length(Win_vec)*length(Enhan_vec);
time_index = zeros(MAXGEN,1);
Obj_index = zeros(MAXGEN,1);
Gen = 0;
for Num_i= 1:1:length(Fea_vec) % searching range for hyperparameter 1

    for Num_j=1:1:length(Win_vec) % searching range for hyperparameter 2

        for Num_k=1:1:length(Enhan_vec) % range for hyperparameter 3
            NumFea = Fea_vec(Num_i);
            NumWin = Win_vec(Num_j);
            NumEnhan = Enhan_vec(Num_k);
            
            rand('state',1)
            for i=1:NumWin
                WeightFea=2*rand(size(train_x,2)+1,NumFea)-1;
               
                WF{i}=WeightFea;
            end                                                         
            
            WeightEnhan=2*rand(NumWin*NumFea+1,NumEnhan)-1;

            [NetoutValidation,test_error_Y] = bls_train_Y_noTest(...
                train_x,train_y,validation_x,WF,WeightEnhan,s,C,NumFea,...
                NumWin,Vp_true,den_true,f,modes_num_vec,index_vec,...
                index_vec_all,mean_y,std_y,validation_x_norm);

            result = [result; NumFea NumWin NumEnhan test_error_Y ...
                NetoutValidation]; % recording all the searching reaults
            
            [Y,I] = min(result(:,4));
            optimal_value = result(I,5:end);
            Obj_index(Gen+1) = Y;
            time_index(Gen+1) = toc;
            Gen = Gen+1;
        end

    end

end
all_index = [time_index Obj_index];
%% optimal prediction
temp_index = find(result(:,4)==min(result(:,4)));
NumFea_hat = result(temp_index,1);
NumWin_hat = result(temp_index,2);
NumEnhan_hat = result(temp_index,3);

Y_hat = result(temp_index,5:end);

end

