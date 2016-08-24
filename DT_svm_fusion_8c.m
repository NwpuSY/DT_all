% function [ Out_label, Out_score] = DT_svm_fusion( K,L,classID,DTL,A_prob,B_prob,C_prob,D_prob,E_prob,F_prob )
function [ Out_label, Out_score] = DT_svm_fusion_8c( M,L,classID,DTL,A_prob,B_prob,C_prob,D_prob,A_prob2,B_prob2,C_prob2,D_prob2 )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
sample_num = size(classID,1);
Out_label =zeros(sample_num,1);
Out_score =zeros(sample_num,1);
cnt=0;

%%д���������ľ������棬���ݴ˽��з����ж�
for i=1:sample_num

   DPinput=[A_prob(i);B_prob(i);C_prob(i);D_prob(i);A_prob2(i);B_prob2(i);C_prob2(i);D_prob2(i)];
   
% % % %    PP = DTL{1};
% % % %    disp_pos = abs(DPinput-PP);
% % % %    de(1) = min(disp_pos);
% % % %    PP=DTL{2};    
% % % %    disp_pos = abs(DPinput-PP);
% % % %    de(2) = min(disp_pos);  
   PP=DTL{1};
   de(1)=1-(sum(sum((DPinput-PP).^2)))/L;
   PP=DTL{2};  
   de(2)=1-(sum(sum((DPinput-PP).^2)))/L;

   [num,val]=sort(de);   %%%�����С�������� 
%   Out_score(i) = 1/(1+de(1));     %2
 	Out_score(i) = de(1)/sum(de);     %1   
   
   if val(1) ==1
       Out_label(i) = 1;
%       Out_score(i) = 1-de(val(1))/sum(de);     %1
%          Out_score(i) = max(DPinput);          %3

        cnt=cnt+1;
   else
       Out_label(i) = 0;
       %     Out_score(i) = 1-de(val(2))/sum(de);   %1
%         Out_score(i) = min(DPinput);          %3

   end
end


% %% ����ָ��


end

