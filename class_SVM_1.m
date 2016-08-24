function [Outputs_,Trn_score_1,cnt,record_zero_line]= class_SVM_1(DTI_Test,ID_TrainDrug,ID_TestDrug,Kd,Kt);
    cutoff_Super = 1.1; % default 1.1
    SuperTargetDTI = GetSuperTargetDTI(DTI_Test,Kt,cutoff_Super);
    useKNN = 1;cnt=0;record_zero_line = [];
    %Output individual
    TrnDTI_S2 = DTI_Test(ID_TrainDrug,:);
    TstDTI_S2 = DTI_Test(ID_TestDrug, :);    % find testing entries of both interactions and non-interactions.

    Trn_score_1 = PredictMonopartite(TrnDTI_S2,TrnDTI_S2, Kd(ID_TrainDrug,ID_TrainDrug), Kd(ID_TrainDrug,ID_TrainDrug), useKNN);
    Out_ = PredictMonopartite(TrnDTI_S2,TstDTI_S2, Kd(ID_TrainDrug,ID_TrainDrug), Kd(ID_TestDrug,ID_TrainDrug), useKNN);
    
   
    %Output super
    TrnDTI_S2_s = SuperTargetDTI(ID_TrainDrug,:);
    TstDTI_S2_s =  SuperTargetDTI (ID_TestDrug , :) ; % find testing entries of both interactions and non-interactions.
    Out_s = PredictMonopartite(TrnDTI_S2_s,TstDTI_S2_s, Kd(ID_TrainDrug,ID_TrainDrug), Kd(ID_TestDrug,ID_TrainDrug),useKNN);
    
    %Final scores and original labels
    Outputs_ =(Out_s .* Out_); % must be values >0
% %     OriginalOutput_{k_fold_r}= TstDTI_S2;
% %     
% %     TrueScore = Outputs_{k_fold_r}(TstDTI_S2==1);
% %     FalseScore= Outputs_{k_fold_r}(TstDTI_S2~=1);
    
   
% %     % Assess the performance in the current CV
% %     [AUC_S2(k_fold_r), AUPR_S2(k_fold_r) ]=EstimationAUC(TrueScore,FalseScore,2000,0);
% %     disp('-');
     for kk=1:size(TrnDTI_S2,2)
        if sum(TrnDTI_S2(:,kk))==0
                cnt =cnt+1;
               record_zero_line =[record_zero_line kk];
        end
    end
end

function [Outputs_, areaROC_ , areaPR_ ] = PredictMonopartite(TrnDTI,TstDTI,TrnSimlarity,TstSimlarity, useKNN)
% averaged assessment from drug or target at one time
if nargin<5
    useKNN = true;
end

if useKNN
    Num=3;Smooth=1; % default parameters and requirements of MLKNN
    TrnDTI = TrnDTI'; %transpose for MLKNN
    TstDTI= TstDTI';%transpose for MLKNN
    
    TrnDTI( TrnDTI~=1)  =-1; %required by MLKNN
    TstDTI( TstDTI~=1)  =-1; %required by MLKNN
    
    %KNN Classifier
    [Prior,PriorN,Cond,CondN]=MLKNN_trainWithKernel(TrnSimlarity,TrnDTI,Num,Smooth);
%     [~,~,~,~,~,Outputs_,~]...
%         =MLKNN_testWithKernel(TrnDTI,TstSimlarity,TstDTI,Num,Prior,PriorN,Cond,CondN);
    [Outputs_,Pre_Labels]...
        =MLKNN_testWithKernel(TrnDTI,TstSimlarity,TstDTI,Num,Prior,PriorN,Cond,CondN);
    
%     Outputs_= Outputs_'; %transpose
    %
%     TstDTI=TstDTI'; %transpose again
    
else
    % We can use other classfiers here. %
    
end
% % % % %6 split scores
% % % % TrueScore = Outputs_(TstDTI==1);
% % % % FalseScore= Outputs_(TstDTI~=1);
% % % % % disp('S2/S3 in S4')
% % % % [areaROC_ , areaPR_ ]=EstimationAUC(TrueScore,FalseScore,2000,0);
end