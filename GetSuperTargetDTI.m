function SuperTargetDTI= GetSuperTargetDTI(DTI,Kt, cutoff)
% T: clustering results containing cluster labels of columns
DTI(DTI ~=1)= 0;

% # --- clustering
dissimilarity = 1- Kt;
Y=squareform(dissimilarity-diag(diag(dissimilarity)) );
Z=linkage(Y,'ward');
if nargin <3
    cutoff = 0.7*(max(Z(:,3)));
end

T = cluster(Z,'cutoff',cutoff,'criterion','distance');

% figure,[H,~] = dendrogram(Z,size(Kt,1),'colorthreshold',cutoff);'default';
% set(H,'LineWidth',2)

% # --- Shrinking
DTI_shrink = DTI;

uT=unique(T);
MergedColIdx  =[];
t=0;
for u=1:length(uT)
    idx = find(T==uT(u));
    if length(idx)<1 %
        continue;
    else
        t=t+1;
        MergedColIdx{t,1} = idx;
    end
end

% # ---Get the interaction matrix between drugs and Super-targets
MergedCol = zeros(size(DTI_shrink,1),length(MergedColIdx));
MergedColCell = cell(1,length(MergedColIdx));

for t=1:length(MergedColIdx)
    colIdx = MergedColIdx{t};
    for c =1: length(colIdx)
        MergedCol(:,t) =  MergedCol(:,t) | DTI_shrink(:, colIdx(c) );
    end
    MergedColCell{t}= repmat(MergedCol(:,t),1,length(colIdx));
end

allMergedCol = cell2mat( MergedColIdx );

%
DTI_shrink(:,allMergedCol) = cell2mat( MergedColCell);

SuperTargetDTI = DTI_shrink;