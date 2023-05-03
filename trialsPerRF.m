% This script is to test the number of trials needed for a flicker
% frequency choice. 

% We will see after how many trials a frequency choice resulted in
% significant (p<0.5) power fold change in a channel in seizure onset zone
% channel. 


%define root dir (where project data is) and sessions of interest:
root_dir='Y:/';
[~,sessions]=fetch_flicker_subjectIDs(root_dir,'flickerneuro');
%[subjectIDs,sessions]=fetch_flicker_subjectIDs(root_dir,'flickerfreq');
%[subjectIDs,sessions]=fetch_flicker_subjectIDs(root_dir,'spep');
%[subjectIDs,sessions]=fetch_flicker_subjectIDs(root_dir,'all');

% TODO: set up a loop to go through all experiment runs
exp_nber=1;

%get soz channels:
fnames=struct;
fnames.root_dir=root_dir;
fnames.subjectID=sessions{exp_nber,'sub'}{:};
fnames.task=sessions{exp_nber,'task'}{:};
fnames.ses=sessions{exp_nber,'ses'}{:};
fetch_channels(fnames,'patho_channels');

%get psd data:
exp_nber=1;
PSD_results=importdata([root_dir '/stg-preproc/sub-' sessions{exp_nber,'sub'}{:} '/task-' sessions{exp_nber,'task'}{:} '/ses-' sessions{exp_nber,'ses'}{:} '/LFP/static_ent/sub-' sessions{exp_nber,'sub'}{:} '_stg-analysis_task-' sessions{exp_nber,'task'}{:} '_ses-' sessions{exp_nber,'ses'}{:} '_nat-psd-refLaplacian.mat']);

%get p-values of fold-change in power for each channel and condition:
p_value_table=readtable([root_dir '/stg-analyses/task-' sessions{exp_nber,'task'}{:} '/sub-' sessions{exp_nber,'sub'}{:} '/ses-' sessions{exp_nber,'ses'}{:} '/LFP/static_ent/LFP_pvalue_table_refLaplacian.csv'],'VariableNamingRule','preserve','RowNamesColumn',1);

% select the channels with p-value<0.05 for 40 Hz AV (2nd column), for 5.5
% Hz is column 5 and 80 Hz is column 8
rows = (p_value_table.(2)<0.05);
p_value_sig_40AV = p_value_table(rows,2);

% check to see after how many trials their p-value became <0.05
% 
