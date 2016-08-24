% function [ P,AAC_prob_animal_new ] = two_label_svm0516( prd_y,predict_label,classID )
function AAC_prob = two_label_svm0516(Trn_Outputs_,train_true_label )

    L=2;
    AAC_prob={};
% % %     ClassID=classID;
   [Xindex,~] = find(train_true_label==1);
   prob = Trn_Outputs_(Xindex,:);
   AAC_prob{1} = prob;
   Xindex= [];
   [Xindex,~] = find(train_true_label==0);
   prob = Trn_Outputs_(Xindex,:);
   AAC_prob{2} = prob;

end

