function [Outputs_ , OriginalOutput, AUC_re, AUPR_re ] = Predicting_DT_main(DTI,SimD_cell, SimT_cell, nCV)

% In the scenario S2 or S3 (dx-ty), generate samples indices for k-fold cross validation
%Example:
%  Predicting_DT_main(DTI,{d_s,d_ATC}, {t_s,t_s_Class}); %for predicting interactions for new drugs,  S2;
%  Predicting_DT_main(DTI', {t_s,t_s_Class},{d_s,d_ATC}); %for predicting interactions for new targets, S3;
%
% INPUT:
% SimD_cell, SimT_cell  -- the cells which store several kinds of
% drug similarities and targets simillarities respectively.
% nCV                           -- the number of fold in cross validation
%
% OUTPUT: 
% OutputMat -- all predicted scores
%
% NOTE:
% This function had been proven that  performing on DTI by PerformS2 is
% equivalent to performing on DTI' by PerformS3. In addition, PerformS3 is
% provided to validate PerformS2 in a symmetric way.
%


if nargin <4
    nCV=5;
end

%make similarity matrices symmetric
for s=1:length(SimD_cell)
    SimD_cell{s}=  (SimD_cell{s}+ SimD_cell{s}' )/2;
end

for s=1:length(SimT_cell)
    SimT_cell{s}=  (SimT_cell{s}+ SimT_cell{s}' )/2;
end


% nCV = 5;
[nDrug, ~]= size(DTI);

% % rand('state',1234567890); % fix random seed for debug, S2
% % CV_Idx_Row= GenerateIdxForCV(nDrug,nCV);
%  2 Predicting

[Outputs_ , OriginalOutput, AUC_re, AUPR_re]= PerformS2(DTI,SimD_cell,SimT_cell,nCV); %
% end


%% Subfunctions
function [Outputs_ , OriginalOutput_, AUC_re, AUPR_re] = PerformS2(DTI,SimD_cell,SimT_cell,nCV)
% Input: CV_Idx_Row= GenerateIdxForCV(nDrug,nCV); %split rows
%     disp('--> supertarget: S2/S3 in S4')
[nDrug, ~]= size(DTI);
rand('state',1234567890); % fix random seed for debug, S2
CV_Idx_Row= GenerateIdxForCV(nDrug,nCV);
cutoff_Super = 1.1; % default 1.1
DTI_Test = DTI;
% % SuperTargetDTI = GetSuperTargetDTI(DTI,Kt,cutoff_Super);

AUC_S2 = zeros(nCV,1);
AUPR_S2 = zeros(nCV,1);
Outputs_= cell(nCV,1); OriginalOutput_= cell(nCV,1);

useSVM = true;
if ~useSVM
    disp('Other Classifier')
end


for k_fold_r = 1:nCV
    cnt=0;record_zero_line = [];
    % find training entries of both interactions and non-interactions
    % Split training set and testing set.
    CV_Temp= CV_Idx_Row;    CV_Temp(k_fold_r) = [];
    ID_TrainDrug =  cell2mat(CV_Temp) ;
    clear CV_Temp;
    ID_TestDrug =  CV_Idx_Row{k_fold_r} ;
    %Output individual
    TrnDTI_S2 = DTI_Test(ID_TrainDrug,:)';
    TstDTI_S2 =  DTI_Test (ID_TestDrug , :)' ; % find testing entries of both interactions and non-interactions.
   
    [Tst_Outputs1_,Trn_Outputs_,cnt1,record_zero_line1]= class_SVM_1(DTI_Test,ID_TrainDrug,ID_TestDrug,SimD_cell{1},SimT_cell{1});

    TrnDTI_S2(record_zero_line1,:)=[];
    TstDTI_S2(record_zero_line1,:)=[];
    Tst_Outputs1_(record_zero_line1,:)=[];
    Trn_Outputs_(record_zero_line1,:)=[];
    
    train_true_label = TrnDTI_S2(:);
    test_true_label = TstDTI_S2(:);
    predict_train_score =  Trn_Outputs_(:);
    predict_test_score1 =  Tst_Outputs1_(:);
    
    A_prob_new1 = two_label_svm0516( predict_train_score,train_true_label );
    
    TrueScore = predict_test_score1(test_true_label==1);
    FalseScore= predict_test_score1(test_true_label~=1);
    
        % Assess the performance in the current CV
    [AUC_S22, AUPR_S22]=EstimationAUC(TrueScore,FalseScore,2000,0);
    disp('-');

    [Tst_Outputs2_,Trn_Outputs_,cnt2,record_zero_line2]= class_SVM_1(DTI_Test,ID_TrainDrug,ID_TestDrug,SimD_cell{1},SimT_cell{2});
    %  change to one line
    Tst_Outputs2_(record_zero_line2,:)=[];
    Trn_Outputs_(record_zero_line2,:)=[];
    
    predict_train_score =  Trn_Outputs_(:);
    predict_test_score2 =  Tst_Outputs2_(:);
    B_prob_new1 = two_label_svm0516( predict_train_score,train_true_label );
    TrueScore = predict_test_score2(test_true_label==1);
    FalseScore= predict_test_score2(test_true_label~=1);
        % Assess the performance in the current CV
    [AUC_S22, AUPR_S22]=EstimationAUC(TrueScore,FalseScore,2000,0);
    disp('-');
    
   
    [Tst_Outputs3_,Trn_Outputs_,cnt3,record_zero_line3]= class_SVM_1(DTI_Test,ID_TrainDrug,ID_TestDrug,SimD_cell{1},SimT_cell{3});
    %  change to one line
     Tst_Outputs3_(record_zero_line3,:)=[];
    Trn_Outputs_(record_zero_line3,:)=[];
    predict_train_score =  Trn_Outputs_(:);
    predict_test_score3 =  Tst_Outputs3_(:);
    C_prob_new1 = two_label_svm0516( predict_train_score,train_true_label );
    TrueScore = predict_test_score3(test_true_label==1);
    FalseScore= predict_test_score3(test_true_label~=1);
        % Assess the performance in the current CV
    [AUC_S22, AUPR_S22]=EstimationAUC(TrueScore,FalseScore,2000,0);
    disp('-');

    [Tst_Outputs4_,Trn_Outputs_,cnt4,record_zero_line4]= class_SVM_1(DTI_Test,ID_TrainDrug,ID_TestDrug,SimD_cell{1},SimT_cell{4});
    %  change to one line
    Tst_Outputs4_(record_zero_line4,:)=[];
    Trn_Outputs_(record_zero_line4,:)=[];
    predict_train_score =  Trn_Outputs_(:);
    predict_test_score4 =  Tst_Outputs4_(:);
    D_prob_new1 = two_label_svm0516( predict_train_score,train_true_label );
    TrueScore = predict_test_score4(test_true_label==1);
    FalseScore= predict_test_score4(test_true_label~=1);
        % Assess the performance in the current CV
    [AUC_S22, AUPR_S22]=EstimationAUC(TrueScore,FalseScore,2000,0);
    disp('-');


    
    [Tst_Outputs12_,Trn_Outputs_,cnt12,record_zero_line12]= class_SVM_1(DTI_Test,ID_TrainDrug,ID_TestDrug,SimD_cell{2},SimT_cell{1});
    %  训练和测试样本的真实标签值
    Tst_Outputs12_(record_zero_line12,:)=[];
    Trn_Outputs_(record_zero_line12,:)=[];
    
    train_true_label = TrnDTI_S2(:);
    test_true_label = TstDTI_S2(:);
    predict_train_score =  Trn_Outputs_(:);
    predict_test_score12 =  Tst_Outputs12_(:);
    
    A_prob_new12 = two_label_svm0516( predict_train_score,train_true_label );
    
    TrueScore = predict_test_score1(test_true_label==1);
    FalseScore= predict_test_score1(test_true_label~=1);
    
        % Assess the performance in the current CV
    [AUC_S22, AUPR_S22]=EstimationAUC(TrueScore,FalseScore,2000,0);
    disp('-');

    
    
    [Tst_Outputs22_,Trn_Outputs_,cnt2,record_zero_line22]= class_SVM_1(DTI_Test,ID_TrainDrug,ID_TestDrug,SimD_cell{2},SimT_cell{2});
    %  change to one line
    Tst_Outputs22_(record_zero_line22,:)=[];
    Trn_Outputs_(record_zero_line22,:)=[];
    
    predict_train_score =  Trn_Outputs_(:);
    predict_test_score22 =  Tst_Outputs22_(:);
    B_prob_new12 = two_label_svm0516( predict_train_score,train_true_label );
    TrueScore = predict_test_score22(test_true_label==1);
    FalseScore= predict_test_score22(test_true_label~=1);
        % Assess the performance in the current CV
    [AUC_S22, AUPR_S22]=EstimationAUC(TrueScore,FalseScore,2000,0);
    disp('-');

    
    [Tst_Outputs32_,Trn_Outputs_,cnt3,record_zero_line32]= class_SVM_1(DTI_Test,ID_TrainDrug,ID_TestDrug,SimD_cell{2},SimT_cell{3});
    %  change to one line
    Tst_Outputs32_(record_zero_line32,:)=[];
    Trn_Outputs_(record_zero_line32,:)=[];
    
    predict_train_score =  Trn_Outputs_(:);
    predict_test_score32 =  Tst_Outputs32_(:);
    C_prob_new12 = two_label_svm0516( predict_train_score,train_true_label );
    TrueScore = predict_test_score32(test_true_label==1);
    FalseScore= predict_test_score32(test_true_label~=1);
        % Assess the performance in the current CV
    [AUC_S22, AUPR_S22]=EstimationAUC(TrueScore,FalseScore,2000,0);
    disp('-');
    
    [Tst_Outputs42_,Trn_Outputs_,cnt4,record_zero_line42]= class_SVM_1(DTI_Test,ID_TrainDrug,ID_TestDrug,SimD_cell{2},SimT_cell{4});
    %  change to one line
    Tst_Outputs42_(record_zero_line42,:)=[];
    Trn_Outputs_(record_zero_line42,:)=[];
    
    predict_train_score =  Trn_Outputs_(:);
    predict_test_score42 =  Tst_Outputs42_(:);
    D_prob_new12 = two_label_svm0516( predict_train_score,train_true_label );
    TrueScore = predict_test_score42(test_true_label==1);
    FalseScore= predict_test_score42(test_true_label~=1);
        % Assess the performance in the current CV
    [AUC_S22, AUPR_S22]=EstimationAUC(TrueScore,FalseScore,2000,0);
    disp('-');
    
    %     save now
    M=2;L=8;  %%%M为样本类别数，即分为多少类,L为分类器数目
    [ DTL ] = DT_svm_8c( M,L,A_prob_new1,B_prob_new1,C_prob_new1,D_prob_new1,A_prob_new12,B_prob_new12,C_prob_new12,D_prob_new12);
    [ Outputs_label, Outputs_score] = DT_svm_fusion_8c( M,L,test_true_label,DTL,predict_test_score1,predict_test_score2,predict_test_score3,...
                                          predict_test_score4,predict_test_score12,predict_test_score22,predict_test_score32,predict_test_score42);
    
    Outputs_{k_fold_r} = Outputs_score;
    OriginalOutput_{k_fold_r}= test_true_label;
    TrueScore = Outputs_score(test_true_label==1);
    FalseScore= Outputs_score(test_true_label~=1);
	
    % Assess the performance in the current CV
    [AUC_S2(k_fold_r), AUPR_S2(k_fold_r) ]=EstimationAUC(TrueScore,FalseScore,2000,0);
    disp('-');
   

    
% % % % ===========================    
end
AUC_re =  mean(AUC_S2);
AUPR_re =  mean(AUPR_S2);
disp([ AUC_re,  AUPR_re ])
    
    
    

%%


