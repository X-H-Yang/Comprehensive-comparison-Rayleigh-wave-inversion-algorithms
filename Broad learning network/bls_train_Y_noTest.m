function [NetoutValidation,test_error_Y] = bls_train_Y_noTest(train_x,...
    train_y,validation_x,WF,WeightEnhan,s,C,NumFea,NumWin,Vp_true,...
    den_true,f,modes_num_vec,index_vec,index_vec_all,mean_y,std_y,...
    validation_x_norm)
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

%% training
% % tic
H1 = [train_x,  0.1 * ones(size(train_x,1),1)];
y=zeros(size(train_x,1),NumWin*NumFea);
for i=1:NumWin
    WeightFea=WF{i};
    A1 = H1 * WeightFea;
    A1 = mapminmax(A1);
    clear WeightFea;
    WeightFeaSparse  = sparse_bls(A1,H1,1e-3,50)';
    WFSparse{i}=WeightFeaSparse;
    
    T1 = H1 * WeightFeaSparse;
    [T1,ps1]  =  mapminmax(T1',0,1);
    T1 = T1';
    
    ps(i)=ps1;
    y(:,NumFea*(i-1)+1:NumFea*i)=T1;
end

clear H1;
clear T1;
H2 = [y,  0.1 * ones(size(y,1),1)];

T2 = H2 * WeightEnhan;


T2 = tansig(T2);
T3=[y T2];
clear H2;
clear T2;

WeightTop = (T3'  *  T3+eye(size(T3',1)) * (C)) \ ( T3'  *  train_y);


NetoutTrain = T3 * WeightTop;
clear T3;

%% validation
test_x = validation_x_norm;
HH1 = [test_x .1 * ones(size(test_x,1),1)];
yy1=zeros(size(test_x,1),NumWin*NumFea);
for i=1:NumWin
    WeightFeaSparse=WFSparse{i};ps1=ps(i);
    TT1 = HH1 * WeightFeaSparse;
    TT1  =  mapminmax('apply',TT1',ps1)';
    
    clear WeightFeaSparse; clear ps1;
    yy1(:,NumFea*(i-1)+1:NumFea*i)=TT1;
end
clear TT1;clear HH1;
HH2 = [yy1 .1 * ones(size(yy1,1),1)];
TT2 = tansig(HH2 * WeightEnhan);
TT3=[yy1 TT2];
clear HH2;clear b2;clear TT2;

NetoutValidation = TT3 * WeightTop;

clear TT3;

%% error calculation
Vs_inverted_profile_norm = NetoutValidation(:)';
% rescale
Vs_inverted_profile = zeros(size(Vs_inverted_profile_norm,1),size(Vs_inverted_profile_norm,2));
for i = 1:1:size(Vs_inverted_profile_norm,2)
    Vs_inverted_profile(:,i) = (Vs_inverted_profile_norm(:,i)*std_y(i)) + mean_y(i);
end
invertedDispersions = calDispersions_2(Vs_inverted_profile,Vp_true,den_true,f,modes_num_vec,index_vec_all);
test_error_Y = sqrt(sum((invertedDispersions-validation_x).^2/length(invertedDispersions)));
NetoutValidation = Vs_inverted_profile;