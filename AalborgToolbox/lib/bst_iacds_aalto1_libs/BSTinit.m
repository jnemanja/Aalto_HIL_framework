clc;
% get path of bst folder 
tmp         = mfilename('fullpath');
idx         = regexp(tmp,'\');
bst_home_folder  = tmp(1:idx(end));

% change to bst folder 
cd(bst_home_folder);

% check subfolder:
if(~exist('BST_BlockSet','dir')) 
    disp('----------------------------------------------------------------');
    disp(' Blockset folder not found!');    
end

bst_blockset_folder = [bst_home_folder 'BST_BlockSet\'];

% add paths
disp('----------------------------------------------------------');
disp('--> Adding BST folders...');
try 
    addpath(bst_blockset_folder);
end
disp('--> OK...');
disp('----------------------------------------------------------');


