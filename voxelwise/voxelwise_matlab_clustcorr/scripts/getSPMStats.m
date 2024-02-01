%% Section 0
projDir = '/nobackup/h_vmac/user/robbwh/voxelwise/fdg_pet/pre_post'
%specify the height threshold 
% for quite strict 2 sided test (0.001 uncorr p, others suggest 0.002 (when using AFNI) fsl default 0.01 (not suggested anymore) 
% ptoz 0.001 -2 (command line fsl) = 3.290527
% qt(pnorm(3.290527),df=66) (in R) = 3.444136 
hThr =  3.444136
nImages=136
dof=66
%% 
%%%Section 1
clc; clearvars -except projDir hThr dof nImages; close all;
addpath(genpath('/nobackup/h_vmac/spm12'));
%%
%Section 2
%SPM estimate smoothness needs filenames in a char array
path = strcat(projDir,'/processingInR/restructured_images');
files = dir(strcat(path,'/vol*.nii'));
for i = 1:length(files)
resFilenames{i} = strcat(path,'/',files(i).name);    
end
%%
%Section 3
currVol = 0;
for i = 1:length(resFilenames)
    nVolsInFile = length(spm_vol(resFilenames{i}));
    for j = 1:nVolsInFile
        inVolFns{currVol+j} = strcat(resFilenames{i},',',num2str(j));
    end
    currVol=currVol+nVolsInFile;
end

inVolFns = char(inVolFns);
%%
%Section 4
%Analysis information
% images in volume
% degrees of freedom of analysis

%%
%Section 5
%Import mask image
maskImgFilename = strcat(projDir, '/processingInR/explicit_mask.nii');
maskImgHdr = spm_vol(maskImgFilename);
maskImg = spm_read_vols(maskImgHdr);
%%
%Section 6
%Get Some image properties:

%Get Smoothness information 
[FWHM,VRpv,R] = spm_est_smoothness(inVolFns,maskImgFilename,[nImages dof]);
%Get voxel to mm matrix transformt
M = maskImgHdr.mat(1:3,1:3);
%Get voxel Dimensions
VOX  = sqrt(diag(M'*M))'; 


%%
%Section 7
%Import T Stat Img
tStatFilename = strcat(projDir, '/processingInR/reconstructed_images/tStatWithNas.nii');
tStatImg = spm_vol(tStatFilename);
%Get matrix of voxel indicies 3xS matrix where S is the number of voxels
[x,y,z]        = ind2sub(maskImgHdr.dim,find(maskImg));
XYZ            = [x y z]';

%%
%Section 8
%Get the voxel values at each voxel index within mask.
%We are using T - statistic images in this example.
%In order to keep the script in similar format to spm_getSPM
%We keep the loaded image in a variable called Z

%Uncomment below for positive association
Z = spm_data_read(tStatImg,'xyz',XYZ);

%Uncomment below for negative association
%Z = -1.*spm_data_read(tStatImg,'xyz',XYZ);

voxCount = length(Z);
%Copy of unmodified XYZ and Z since these matricies will undergo height and
%extent thresholding
XYZum = XYZ;
Zum = Z;

%%
%Section 9
%Find which voxels pass a height threshold. Defining this in terms of T
%Default in fsl is a Z-statistic=2.3, default in SPM is a p-value=0.001
indPassingHThr = find(Z > hThr);
Z = Z(:,indPassingHThr);
XYZ = XYZ(:,indPassingHThr);

%%
%Section 10
%Split remaining voxels into clusters.

% FORMAT [N Z M A XYZ] = spm_max(X,L)
% inputs:
% Z     - values of 3-D field
% XYZ     - locations [x y z]' {in voxels}
%
% outputs:
% N     - size of region {in voxels)
% Z     - Z values of maxima
% XYZ     - location of maxima {in voxels}
% A     - region number
% L   - cell array of voxel locations for each cluster
[N,Z,XYZ,A,L]  = spm_max(Z,XYZ);
%Copy cluster values and XYZ
Zexcursion = Z;
XYZexcursion = XYZ;
%number of clusters
c=max(A);

%%
%Section 11
%cluster extent in resels
FWHM = full(FWHM);
V2R   = 1/prod(FWHM);
K = N*V2R;
%search volume in voxels
S = voxCount;
%Script is setup for T-statistic
STAT='T';
n_conj = 1;
PkVec = [];
PnVec = [];
QcVec = [];
QpVec = [];
clustSize = [];
maxStat = [];
maxCoord = [];
locOfVoxInCluster = {};

%%
%Section 12
 while numel(find(isfinite(Z)))
    %-Find largest remaining local maximum
    %------------------------------------------------------------------
    [U,i]  = max(Z);            %-largest maxima
    j      = find(A == A(i));   %-maxima in cluster
    %N(find(A==i)) % Cluster Size 

    %-Compute cluster {k} and peak-level {u} p-values for this cluster
    %------------------------------------------------------------------
    if STAT ~= 'P'

        % p-values
        %--------------------------------------------------------------
        
        % p value of finding 1 or more clusters
        %with a size > 0, given height threshold hThr and assuming the
        %image is a single resel
        %
        %p value assuming cluster act as a poisson point process
        %Pz      = spm_P(1,0,   hThr,[1 dof],STAT,1,n_conj,S);  % uncorrected p value
        %
        % p value of finding 1 or more clusters with a size > 0,
        %given height threshold hThr, corrected for the number of estimated
        %resels
        %
        %p value assuming cluster act as a poisson point process
        %Pu      = spm_P(1,0,   hThr,[1 dof],STAT,R,n_conj,S);
        
        
        % p value of finding 1 or more clusters with a size > K(i) in
        % resels, given height threshold hThr, and 
        % corrected for the number of estimated resels in volume.
        %
        % Pk is the p value FWE-corrected by assuming clusters act
        % as a poisson point process
        %
        %Pn is the raw p value given the criteria above
        %approximated by finding the expected Euler Characteristic as
        %described here: https://www.fil.ion.ucl.ac.uk/spm/course/slides02/overview/Stats.htm 
        [Pk,Pn] = spm_P(1,K(i),hThr,[1 dof],STAT,R,n_conj,S);
        PkVec = [PkVec Pk];
        PnVec = [PnVec Pn];
        clustSize = [clustSize N(i)];
        maxStat = [maxStat U];
        maxCoord = [maxCoord XYZ(:,i)];
        locOfVoxInCluster{end + 1} = L{A(i)};
    end
    % Set local maxima to NaN
    Z(j) = NaN;
    
 end
            
%%
%Section 13
% q-values (FDR)
%--------------------------------------------------------------

[~,ind]=sort(PnVec);

PclustFDR = spm_P_FDR(sort(PnVec));

unsorted = 1:length(PnVec);

revInd(ind) = unsorted;

OrigOrderPclustFDR = PclustFDR(revInd);

%%
%Section 14
%Reorder statistics by cluster size

[sortedClustSize, idx] = sort(clustSize,'descend');

sortedClusterIdx = 1:length(sortedClustSize);

sortedPkVec = PkVec(idx);

sortedPnVec = PnVec(idx);

sortedPclustFDR = OrigOrderPclustFDR(idx);

sortedMaxCoord = maxCoord(:,idx);

lenMaxCoord = length(sortedMaxCoord)

if lenMaxCoord ==0
	disp('hurr')
else



%Subtract 1 to match fsleyes viewing from 0 to n-1 in a given dimension
sortedXCoord = sortedMaxCoord(1,:) - 1;
sortedYCoord = sortedMaxCoord(2,:) - 1;
sortedZCoord = sortedMaxCoord(3,:) - 1;

sortedMaxStat = maxStat(idx);
sortedLocOfVoxInCluster = locOfVoxInCluster(idx);

statsTableColumnNames = {'ClusterNumber','ClusterSizeInVoxels','ClusterFWEcorrp','ClusterUncorrp','ClusterFDRcorrp','PeakTstatistic','PeakVoxelCoorX', 'PeakVoxCoorY','PeakVoxCoorZ'};

statsTable = table(sortedClusterIdx(:),sortedClustSize(:),sortedPkVec(:),sortedPnVec(:),sortedPclustFDR(:),sortedMaxStat(:),sortedXCoord(:),sortedYCoord(:),sortedZCoord(:));
statsTable.Properties.VariableNames = statsTableColumnNames;

%%
%Section 15
statsTableFn = strcat(projDir, '/processingInR/spm_results/myposResults.txt');
writetable(statsTable,statsTableFn,'Delimiter',' ')  
end
%%
%Section 16
%Vector of cluster numbers you want an outputted mask from the statsTable
%specify clustersToOutput as a matlab vector:
%Example output clusters 1,5,7: clustersToOutput = [1 5 7];
%Example output clusters 1,2,3,4,5: clustersToOutput = 1:5;

if lenMaxCoord ==0
	disp('hurr')
else

clustersToOutput = 1:1;
outpath=strcat(projDir, '/processingInR/spm_results');
tStatsImgData = spm_read_vols(tStatImg);
if length(clustersToOutput) > 0
    for i = 1:length(clustersToOutput)
        clusterMaskFn = strcat(outpath,'/my_pos_cluster_',num2str(clustersToOutput(i)),'_posmask.nii')
        clusterMask = zeros(size(tStatsImgData));
        clusterMask(sub2ind(size(tStatsImgData), sortedLocOfVoxInCluster{clustersToOutput(i)}(1,:),sortedLocOfVoxInCluster{clustersToOutput(i)}(2,:),sortedLocOfVoxInCluster{clustersToOutput(i)}(3,:))) = 1;
        clusterMaskHdr = maskImgHdr;
        clusterMaskHdr.fname = clusterMaskFn;
        clusterMaskHdr.descrip = 'cluster mask built in matlab using spm12';
        clusterMaskHdr = rmfield(clusterMaskHdr,'private');
        clusterMaskHdr = spm_write_vol(clusterMaskHdr,clusterMask);
    end
end
end
%% 
%%%Section 1
clc; clearvars -except projDir hThr dof nImages; close all;
addpath(genpath('/nobackup/h_vmac/spm12'));
%%
%Section 2
%SPM estimate smoothness needs filenames in a char array
path = strcat(projDir,'/processingInR/restructured_images');
files = dir(strcat(path,'/vol*.nii'));
for i = 1:length(files)
resFilenames{i} = strcat(path,'/',files(i).name);    
end
%%
%Section 3
currVol = 0;
for i = 1:length(resFilenames)
    nVolsInFile = length(spm_vol(resFilenames{i}));
    for j = 1:nVolsInFile
        inVolFns{currVol+j} = strcat(resFilenames{i},',',num2str(j));
    end
    currVol=currVol+nVolsInFile;
end

inVolFns = char(inVolFns);
%%
%Section 4
%Analysis information
% images in volume
% degrees of freedom of analysis

%%
%Section 5
%Import mask image
maskImgFilename = strcat(projDir, '/processingInR/explicit_mask.nii');
maskImgHdr = spm_vol(maskImgFilename);
maskImg = spm_read_vols(maskImgHdr);
%%
%Section 6
%Get Some image properties:

%Get Smoothness information 
[FWHM,VRpv,R] = spm_est_smoothness(inVolFns,maskImgFilename,[nImages dof]);
%Get voxel to mm matrix transformt
M = maskImgHdr.mat(1:3,1:3);
%Get voxel Dimensions
VOX  = sqrt(diag(M'*M))'; 


%%
%Section 7
%Import T Stat Img
tStatFilename = strcat(projDir, '/processingInR/reconstructed_images/tStatWithNas.nii');
tStatImg = spm_vol(tStatFilename);
%Get matrix of voxel indicies 3xS matrix where S is the number of voxels
[x,y,z]        = ind2sub(maskImgHdr.dim,find(maskImg));
XYZ            = [x y z]';

%%
%Section 8
%Get the voxel values at each voxel index within mask.
%We are using T - statistic images in this example.
%In order to keep the script in similar format to spm_getSPM
%We keep the loaded image in a variable called Z

%Uncomment below for positive association
%Z = spm_data_read(tStatImg,'xyz',XYZ);

%Uncomment below for negative association
Z = -1.*spm_data_read(tStatImg,'xyz',XYZ);

voxCount = length(Z);
%Copy of unmodified XYZ and Z since these matricies will undergo height and
%extent thresholding
XYZum = XYZ;
Zum = Z;

%%
%Section 9
%Find which voxels pass a height threshold. Defining this in terms of T
%Default in fsl is a Z-statistic=2.3, default in SPM is a p-value=0.001
indPassingHThr = find(Z > hThr);
Z = Z(:,indPassingHThr);
XYZ = XYZ(:,indPassingHThr);

%%
%Section 10
%Split remaining voxels into clusters.

% FORMAT [N Z M A XYZ] = spm_max(X,L)
% inputs:
% Z     - values of 3-D field
% XYZ     - locations [x y z]' {in voxels}
%
% outputs:
% N     - size of region {in voxels)
% Z     - Z values of maxima
% XYZ     - location of maxima {in voxels}
% A     - region number
% L   - cell array of voxel locations for each cluster
[N,Z,XYZ,A,L]  = spm_max(Z,XYZ);
%Copy cluster values and XYZ
Zexcursion = Z;
XYZexcursion = XYZ;
%number of clusters
c=max(A);

%%
%Section 11
%cluster extent in resels
FWHM = full(FWHM);
V2R   = 1/prod(FWHM);
K = N*V2R;
%search volume in voxels
S = voxCount;
%Script is setup for T-statistic
STAT='T';
n_conj = 1;
PkVec = [];
PnVec = [];
QcVec = [];
QpVec = [];
clustSize = [];
maxStat = [];
maxCoord = [];
locOfVoxInCluster = {};

%%
%Section 12
 while numel(find(isfinite(Z)))
    %-Find largest remaining local maximum
    %------------------------------------------------------------------
    [U,i]  = max(Z);            %-largest maxima
    j      = find(A == A(i));   %-maxima in cluster
    %N(find(A==i)) % Cluster Size 

    %-Compute cluster {k} and peak-level {u} p-values for this cluster
    %------------------------------------------------------------------
    if STAT ~= 'P'

        % p-values
        %--------------------------------------------------------------
        
        % p value of finding 1 or more clusters
        %with a size > 0, given height threshold hThr and assuming the
        %image is a single resel
        %
        %p value assuming cluster act as a poisson point process
        %Pz      = spm_P(1,0,   hThr,[1 dof],STAT,1,n_conj,S);  % uncorrected p value
        %
        % p value of finding 1 or more clusters with a size > 0,
        %given height threshold hThr, corrected for the number of estimated
        %resels
        %
        %p value assuming cluster act as a poisson point process
        %Pu      = spm_P(1,0,   hThr,[1 dof],STAT,R,n_conj,S);
        
        
        % p value of finding 1 or more clusters with a size > K(i) in
        % resels, given height threshold hThr, and 
        % corrected for the number of estimated resels in volume.
        %
        % Pk is the p value FWE-corrected by assuming clusters act
        % as a poisson point process
        %
        %Pn is the raw p value given the criteria above
        %approximated by finding the expected Euler Characteristic as
        %described here: https://www.fil.ion.ucl.ac.uk/spm/course/slides02/overview/Stats.htm 
        [Pk,Pn] = spm_P(1,K(i),hThr,[1 dof],STAT,R,n_conj,S);
        PkVec = [PkVec Pk];
        PnVec = [PnVec Pn];
        clustSize = [clustSize N(i)];
        maxStat = [maxStat U];
        maxCoord = [maxCoord XYZ(:,i)];
        locOfVoxInCluster{end + 1} = L{A(i)};
    end
    % Set local maxima to NaN
    Z(j) = NaN;
    
 end
            
%%
%Section 13
% q-values (FDR)
%--------------------------------------------------------------

[~,ind]=sort(PnVec);

PclustFDR = spm_P_FDR(sort(PnVec));

unsorted = 1:length(PnVec);

revInd(ind) = unsorted;

OrigOrderPclustFDR = PclustFDR(revInd);

%%
%Section 14
%Reorder statistics by cluster size

[sortedClustSize, idx] = sort(clustSize,'descend');

sortedClusterIdx = 1:length(sortedClustSize);

sortedPkVec = PkVec(idx);

sortedPnVec = PnVec(idx);

sortedPclustFDR = OrigOrderPclustFDR(idx);

sortedMaxCoord = maxCoord(:,idx);
%Subtract 1 to match fsleyes viewing from 0 to n-1 in a given dimension
sortedXCoord = sortedMaxCoord(1,:) - 1;
sortedYCoord = sortedMaxCoord(2,:) - 1;
sortedZCoord = sortedMaxCoord(3,:) - 1;

sortedMaxStat = maxStat(idx);
sortedLocOfVoxInCluster = locOfVoxInCluster(idx);

statsTableColumnNames = {'ClusterNumber','ClusterSizeInVoxels','ClusterFWEcorrp','ClusterUncorrp','ClusterFDRcorrp','PeakTstatistic','PeakVoxelCoorX', 'PeakVoxCoorY','PeakVoxCoorZ'};

statsTable = table(sortedClusterIdx(:),sortedClustSize(:),sortedPkVec(:),sortedPnVec(:),sortedPclustFDR(:),sortedMaxStat(:),sortedXCoord(:),sortedYCoord(:),sortedZCoord(:));
statsTable.Properties.VariableNames = statsTableColumnNames;

%%
%Section 15
statsTableFn = strcat(projDir, '/processingInR/spm_results/mynegResults.txt');
writetable(statsTable,statsTableFn,'Delimiter',' ')  
%%
%Section 16
%Vector of cluster numbers you want an outputted mask from the statsTable
%specify clustersToOutput as a matlab vector:
%Example output clusters 1,5,7: clustersToOutput = [1 5 7];
%Example output clusters 1,2,3,4,5: clustersToOutput = 1:5;
clustersToOutput = 1:5;
outpath=strcat(projDir, '/processingInR/spm_results');
tStatsImgData = spm_read_vols(tStatImg);
if length(clustersToOutput) > 0
    for i = 1:length(clustersToOutput)
        clusterMaskFn = strcat(outpath,'/my_neg_cluster_',num2str(clustersToOutput(i)),'_negmask.nii')
        clusterMask = zeros(size(tStatsImgData));
        clusterMask(sub2ind(size(tStatsImgData), sortedLocOfVoxInCluster{clustersToOutput(i)}(1,:),sortedLocOfVoxInCluster{clustersToOutput(i)}(2,:),sortedLocOfVoxInCluster{clustersToOutput(i)}(3,:))) = 1;
        clusterMaskHdr = maskImgHdr;
        clusterMaskHdr.fname = clusterMaskFn;
        clusterMaskHdr.descrip = 'cluster mask built in matlab using spm12';
        clusterMaskHdr = rmfield(clusterMaskHdr,'private');
        clusterMaskHdr = spm_write_vol(clusterMaskHdr,clusterMask);
    end
end
