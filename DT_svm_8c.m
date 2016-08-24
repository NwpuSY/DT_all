% function [ DTL ] = DT_svm( L,prob_animal_new1,prob_animal_new2,prob_animal_new3,prob_animal_new4,prob_animal_new5,prob_animal_new6)
function [ DTL ] = DT_svm_8c( M,L,prob_animal_new1,prob_animal_new2,prob_animal_new3,prob_animal_new4,prob_animal_new12,prob_animal_new22,prob_animal_new32,prob_animal_new42)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% load('GreatID_DeleteRepeat.mat')

DTL={};
for i=1:M
    AAC_1 = prob_animal_new1{i};
    AAC_2 = prob_animal_new2{i};
    AAC_3 = prob_animal_new3{i};
    AAC_4 = prob_animal_new4{i};
    AAC_12 = prob_animal_new12{i};
    AAC_22 = prob_animal_new22{i};
    AAC_32 = prob_animal_new32{i};
    AAC_42 = prob_animal_new42{i};
    
    
    
%     PSSM_DP = PSSM_prob_animal_new{i};
    len=length(AAC_1);
    DP_k=cell(len,1);
    for j=1:len
        DP_k{j}=[AAC_1(j,:);AAC_2(j,:);AAC_3(j,:);AAC_4(j,:);AAC_12(j,:);AAC_22(j,:);AAC_32(j,:);AAC_42(j,:)];
    end
    DT1=zeros(L,1);
    for j=1:len
        PP=DP_k{j};
        DT1 = DT1+PP;
    end   
    DT1=DT1/len;
    DTL{i}=DT1;
end




end

