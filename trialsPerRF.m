% This script is to test the number of trials needed for a flicker
% frequency choice. 

% We will see after how many trials a frequency choice resulted in
% significant (p<0.5) power fold change in a channel in seizure onset zone
% channel. 


%define root dir (where project data is) and sessions of interest:
root_dir='Y:\';
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
fnames.analysis_folder = 'stg-analysis';

patho_channels=fetch_channels(fnames,'patho_channels');
contain = 'soz';
row_path = (patho_channels{:,'feature'});
soz_rows = contains(row_path, contain);
soz_channels=patho_channels(soz_rows,["label","feature"]);

%get psd data: TODO - can we just get the soz channels psd?
% filter based on labels (channels), such that are in soz.
exp_nber=1;
data = [root_dir 'stg-preproc\sub-' sessions{exp_nber,'sub'}{:}...
    '\task-' sessions{exp_nber,'task'}{:} '\ses-' sessions{exp_nber,'ses'}{:}...
    '\LFP\static_ent\sub-' sessions{exp_nber,'sub'}{:} '_stg-analysis_task-'...
    sessions{exp_nber,'task'}{:} '_ses-' sessions{exp_nber,'ses'}{:}...
    '_nat-psd-refLaplacian.mat'];
PSD_results=importdata(data);

% get the soz channels PSD
% selected_rows=cell(size(soz_channels{:,"label"},1),1);
% psd_results_soz_ch=cell(size(soz_channels{:,"label"},1),1);
% 
% for str=1:size(soz_channels{:,"label"},1)
%     selected_rows =  contains(PSD_results.label,soz_channels{str,"label"});
%     psd_results_soz_ch(str) = {PSD_results.label(selected_rows)};
% end
% psd_results_soz_ch_vertcat = vertcat(psd_results_soz_ch{:});
% psd_results_soz_ch_mat = unique(psd_results_soz_ch_vertcat);
selected_rows =  contains(PSD_results.label,soz_channels{:,"label"});
psd_results_soz_ch = PSD_results.label(selected_rows);

%get p-values of fold-change in power for each channel and condition: TODO
%- can we calculate p-value per trial
p_value_table=readtable([root_dir 'stg-analyses\task-' sessions{exp_nber,'task'}{:}...
    '\sub-' sessions{exp_nber,'sub'}{:} '\ses-' sessions{exp_nber,'ses'}{:}...
    '\LFP\static_ent\LFP_pvalue_table_refLaplacian.csv'],'VariableNamingRule',...
    'preserve','RowNamesColumn',1);

% select the channels with p-value<0.05 for 40 Hz AV (2nd column), for 5.5
% Hz is column 5 and 80 Hz is column 8
% TODO: can we see after how many trial each freq had a p<0.05?
rows = (p_value_table.(2)<0.05);
p_value_sig_40AV = p_value_table(rows,2);
p_value_rows=matches(p_value_sig_40AV.Row,psd_results_soz_ch,'IgnoreCase',true);
p_value_sig_40AV_soz = p_value_sig_40AV(p_value_rows,:);
p_value_sig_40AV_soz_chs = p_value_sig_40AV_soz.Row(:);

% find the index of PSD_results.labels that corespond to the chs with
% significant p-values in soz


selected_PSD_result_chs = matches(PSD_results.label,p_value_sig_40AV_soz_chs(:));
PSD_results_label_sig_soz_chs = PSD_results.label(selected_PSD_result_chs);
idx_PSD_results_label_sig_soz_ch = find(selected_PSD_result_chs);

% check to see after how many trials their p-value became <0.05
% 
% evaluate_ent_degree_relpower_perTrial(fnames);
% conditions_of_interest=sort(PSD_results.condition(contains(PSD_results.condition,'Hz') & ~contains(PSD_results.condition,'occluded') & ~contains(PSD_results.condition,'min')));
control_condition='Baseline';
conditions_of_interest = {'40Hz-AV'};
num_trials=size(PSD_results.data{1,10}{1,1},1);
% zscore_table=zeros(length(PSD_results.label),length(conditions_of_interest));
% pvalue_table=zeros(length(PSD_results.label),length(conditions_of_interest));
pvalue_trial=zeros(size(p_value_sig_40AV_soz_chs,1),num_trials); %for each channel=zeros(length(PSD_results.label),length(num_trials));
% for i=1:size(zscore_table,2) %for each condition
% 
%     freq_interest=str2num(regexprep(conditions_of_interest{i},'Hz.+',''));
%     [~,temp]=min(abs(PSD_results.data{1,strcmp(PSD_results.condition,conditions_of_interest{i})}{3}-freq_interest)); %find index of frequency closest to frequency of interest- have to do this because sample rate of EDF file not an integer sometimes (error with Natus)
%     freq_interest_index=temp;
i =1;   
for ch = idx_PSD_results_label_sig_soz_ch(1):idx_PSD_results_label_sig_soz_ch(end)
    disp(ch)
   
    stim_values=[];
    baseline_values=[];
    
    for tr=1:num_trials %for however many number of trials of given condition there are
        current_stim_value=PSD_results.data{ch,10}{1,1}(tr,:);
        current_baseline_value=PSD_results.data{ch,2}{1,1}(tr,:);

        stim_values=[stim_values current_stim_value];
        baseline_values=[baseline_values current_baseline_value];
        pvalue_trial(i,tr)=pval_randomshuffle([stim_values' baseline_values'],500);
    end
    i = i + 1;
%     zscore_table(ch,1)=(mean(stim_values)/mean(baseline_values))-1;
%     pvalue_table(ch,1)=pval_randomshuffle([stim_values' baseline_values'],500);
end

% end

%         zscore_table=array2table(zscore_table,'RowNames',PSD_results.label,'VariableNames',conditions_of_interest); %make matrix into table
%         writetable(zscore_table,[fnames.analysis_folder,'/LFP/static_ent/LFP_zscore_table_ref' ref_method{:} '.csv'],'WriteRowNames',1);
%         
%         pvalue_table=array2table(pvalue_table,'RowNames',PSD_results.label,'VariableNames',conditions_of_interest); %make matrix into table
%         writetable(pvalue_table,[fnames.analysis_folder,'/LFP/static_ent/LFP_pvalue_table_ref' ref_method{:} '.csv'],'WriteRowNames',1);
%    
trial(1:num_trials)="trial ";
num = string(1:num_trials);
trial_names=append(trial,num);
pvalue_trial_table=array2table(pvalue_trial,'RowNames',p_value_sig_40AV_soz_chs(:),'VariableNames', trial_names); %make matrix into table
p_value_sig_40AV_soz_chs_trials = pvalue_trial_table(p_value_sig_40AV_soz_chs,:);


% writetable(pvalue_trial_table,[fnames.analysis_folder,'/LFP/static_ent/LFP_pvalue_trial_table_refLaplacian.csv'],'WriteRowNames',1);
       