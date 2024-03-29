% This script is to test the number of trials needed for a flicker
% frequency choice. 

% We will see after how many trials a frequency choice resulted in
% significant (p<0.5) power fold change in a channel in seizure onset zone
% channel. 
clear
close all 
%% set up where the dataset is and folder structure for sessions

%define root dir (where project data is) and sessions of interest:
root_dir='Y:\';
[~,sessions]=fetch_flicker_subjectIDs(root_dir,'flickerneuro');
p_values.ses = join([sessions.sub,sessions.ses],'_',2);

for exp_nber=1:size(p_values.ses,1)-1
    disp("Processing session data: " + sessions{exp_nber,'sub'}{:} + '_ses-'...
        + sessions{exp_nber,'ses'}{:});
    %get soz channels:
    fnames=struct;
    fnames.root_dir=root_dir;
    fnames.subjectID=sessions{exp_nber,'sub'}{:};
    fnames.task=sessions{exp_nber,'task'}{:};
    fnames.ses=sessions{exp_nber,'ses'}{:};
    fnames.analysis_folder = 'stg-analysis';
    
    %% Import the LFP PSD Laplacian refernced data for the channels

    data = [root_dir 'stg-preproc\sub-' sessions{exp_nber,'sub'}{:}...
        '\task-' sessions{exp_nber,'task'}{:} '\ses-' sessions{exp_nber,'ses'}{:}...
        '\LFP\static_ent\sub-' sessions{exp_nber,'sub'}{:} '_stg-analysis_task-'...
        sessions{exp_nber,'task'}{:} '_ses-' sessions{exp_nber,'ses'}{:}...
        '_nat-psd-refLaplacian.mat'];
%     PSD_results=importdata(data);
    PSD_resultsObj= matfile(data); 
    PSD_results = PSD_resultsObj.PSD_results_ref_preproc;   
    %% Get SOZ channels with p-values < 0.05
    
    patho_channels=fetch_channels(fnames,'patho_channels');
    contain = 'soz';
    row_path = (patho_channels{:,'feature'});
    soz_rows = contains(row_path, contain);
    soz_channels=patho_channels(soz_rows,["label","feature"]);
    
    % get the soz channels PSD
    selected_rows =  contains(PSD_results.label,soz_channels{:,"label"});
    psd_results_soz_ch = PSD_results.label(selected_rows);
    
    %get p-values of fold-change in power for each channel and condition: TODO
    %- can we calculate p-value per trial
    p_value_table=readtable([root_dir 'stg-analyses\task-' sessions{exp_nber,'task'}{:}...
        '\sub-' sessions{exp_nber,'sub'}{:} '\ses-' sessions{exp_nber,'ses'}{:}...
        '\LFP\static_ent\LFP_pvalue_table_refLaplacian.csv'],'VariableNamingRule',...
        'preserve','RowNamesColumn',1);
    
    
    %% Nested loops to go through session's conditions and for each channel calculate the p-value for each trial
    
    % check to see after how many trials their p-value became <0.05
    % 
    % evaluate_ent_degree_relpower_perTrial(fnames);
    conditions_of_interest=sort(PSD_results.condition(contains(PSD_results.condition,'Hz-AV') & ~contains(PSD_results.condition,'occluded') & ~contains(PSD_results.condition,'min')));
    p_values.conditions= {};
    control_condition='Baseline';
%     conditions_of_interest = {'80Hz-AV'};
    

    p_values.channels = {};
    p_values.channels.labels = {};
    p_values.channels.run = {};
    p_values.channels.means =[];
    p_values.channels.stdDev =[];
    p_values_comp = table();

    for con=1:size(conditions_of_interest,2) %for each condition
        disp("Processing condition: " + conditions_of_interest{con});
        p_values.conditions(con) = conditions_of_interest(con);
        freq_interest_boolean=strcmp(PSD_results.condition,conditions_of_interest{con}); %find index of frequency closest to frequency of interest- have to do this because sample rate of EDF file not an integer sometimes (error with Natus)
        freq_interest_index=find(freq_interest_boolean);
        disp("freq_interest_index "+ freq_interest_index);
        num_trials=size(PSD_results.data{1,freq_interest_index}{1,1},1);
        baseline_boolean = strcmp(PSD_results.condition,control_condition);
        baseline_index = find(baseline_boolean);
        disp("baseline_condition_index " + baseline_index);
        % select the channels with p-value<0.05 for 40 Hz AV (2nd column), for 5.5
        % Hz is column 5 and 80 Hz is column 8
       
        freq_interest_pvalue_idx = find(strcmp(p_value_table.Properties.VariableNames,conditions_of_interest{con}));
        rows = (p_value_table.(freq_interest_pvalue_idx)<0.05); % after how many trial each freq had a p<0.05
        p_value_sig_condition_of_interest = p_value_table(rows,freq_interest_pvalue_idx);
        p_value_rows=matches(p_value_sig_condition_of_interest.Row,psd_results_soz_ch,'IgnoreCase',true);
        p_value_sig_condition_of_interest_soz = p_value_sig_condition_of_interest(p_value_rows,:);
        p_value_sig_condition_of_interest_soz_chs = p_value_sig_condition_of_interest_soz.Row(:);
        
    
        % find the index of PSD_results.labels that corespond to the chs with
        % significant p-values in soz
        selected_PSD_result_chs = matches(PSD_results.label,p_value_sig_condition_of_interest_soz_chs(:));
        PSD_results_label_sig_soz_chs = PSD_results.label(selected_PSD_result_chs);
        idx_PSD_results_label_sig_soz_ch = find(selected_PSD_result_chs);
        channels(1: size(PSD_results_label_sig_soz_chs,1),1) = {p_values.conditions(con)};
        freq_interest=str2double(regexprep(conditions_of_interest{con},'Hz.+',''));
        [~,temp]=min(abs(PSD_results.data{1,strcmp(PSD_results.condition,conditions_of_interest{con})}{3}-freq_interest)); %find index of frequency closest to frequency of interest- have to do this because sample rate of EDF file not an integer sometimes (error with Natus)
        freq_interest_index_tr=temp;
        disp("The freq index for the trial in session is: "+freq_interest_index_tr)
        trial="";
        trial_names="";
        for ch = 1:size(PSD_results_label_sig_soz_chs,1)
            disp("Processing channel: " + PSD_results_label_sig_soz_chs{ch});
            disp("The channel index in PSD_results is row number: " + idx_PSD_results_label_sig_soz_ch(ch));
            p_values.channels.labels{end+1,1} = {strjoin([channels{ch},PSD_results_label_sig_soz_chs{ch}],'_')};

            for iteration=1:10
                stim_values=[];
                baseline_values=[];
                pvalue_trial=zeros(size(p_value_sig_condition_of_interest_soz_chs,1),num_trials); 
%                 disp("Randomized trial run # " + iteration);
                trial_order= randperm(num_trials);
%                 trial_order = 1:num_trials;
                for tr=trial_order %for however many number of trials of given condition there are
%                     disp("trial number "+tr);
                    current_stim_value=PSD_results.data{idx_PSD_results_label_sig_soz_ch(ch),freq_interest_index}{1,1}(tr,freq_interest_index_tr);
                    current_baseline_value=PSD_results.data{idx_PSD_results_label_sig_soz_ch(ch),baseline_index}{1,1}(tr,freq_interest_index_tr);
            
                    stim_values=[stim_values current_stim_value];
                    baseline_values=[baseline_values current_baseline_value];
                    pvalue_trial(ch,tr)=pval_randomshuffle([stim_values' baseline_values'],1000);
                    
                end
                p_values.channels.run{ch,iteration}=pvalue_trial(ch,:);
            end
          
            % TODO: calculate Mean and std dev per channel
            [l, w ] = size(p_values.channels.run(ch,:));
            pvalues_trial_mat = cell2mat(p_values.channels.run(ch,:));
            pvalues_trial_mat_reshape = transpose(reshape(pvalues_trial_mat,[num_trials, w]));            
            pvalue_trial_mean = mean(pvalues_trial_mat_reshape,1);
            pvalue_trial_stdDev = std(pvalues_trial_mat_reshape,1);
            p_values.channels.means(end+1,1:num_trials)= pvalue_trial_mean;
            p_values.channels.stdDev(end+1,1:num_trials) = pvalue_trial_stdDev;

%             % compare Lou and Sina's trial 15
%             p_values_comp(end+1,:)= array2table([p_values.channels.run{ch,1}(1,num_trials),...
%                 p_value_sig_condition_of_interest{PSD_results_label_sig_soz_chs{ch},1}],...
%                 'VariableNames', ["Sina"; "Lou"]');
%             
%             % PSD plots condition vs. baseline
%             figure("Name",p_values.channels.labels{end,1}{:} )
%             plot_PSD(conditions_of_interest{con},PSD_results_label_sig_soz_chs{ch},PSD_results,PSD_results.label,PSD_results.condition,condition_color('80Hz-AV'),0,1,1)% the 0 is for std dev 
%             hold on;
%             plot_PSD('Baseline',PSD_results_label_sig_soz_chs{ch},PSD_results,PSD_results.label,PSD_results.condition,[0 0 0],0,1,1)
%             xlabel("freq. (Hz)");
%             ylabel("mean power (log_{10})");
%             legend(["stimulation", "baseline"]);
%             hold off;

        end
        
        trial(1:num_trials)="trial ";
        num = string(1:num_trials);
        trial_names=append(trial,num);
        disp("----------");
    end
    %% plotting

    figure("Name",p_values.ses{exp_nber})
    hold on
    semilogy(1:num_trials, p_values.channels.means');
%     plot(1:num_trials, p_values.channels.means');
    errorbar(p_values.channels.means',p_values.channels.stdDev')
    y_p = ones(1,num_trials)*0.05;
    plot(1:num_trials,y_p,"LineStyle","--","LineWidth",2);   
    leg = string(p_values.channels.labels(1:end,1));
    leg_edited = replace(leg,'_','.');
    legend(leg_edited,'Location','southoutside','FontSize',6,'NumColumns',2);
    title_updated = replace(p_values.ses{exp_nber},'_','.');
    title(title_updated)
    xticks(1:1:15)
    set(gca, 'YScale', 'log')
    xlabel("number of trials")
    ylabel("p-value (log_{10})")
    grid on
    hold off
%% Save p_values as .mat file

    p_values_table = array2table(p_values.channels.means ,'RowNames',leg, 'VariableNames', trial_names');%"VariableNames", {colNames},
    dir_path = [root_dir 'Sina\stg-analysis\num_trial_per_freq_choice\'];
    file_path = fullfile(dir_path, [sessions{exp_nber,'sub'}{:} '_ses-' sessions{exp_nber,'ses'}{:}...
        '_LFP_pvalue_trial_table_refLaplacian.mat']);
    fig_path = fullfile(dir_path, [sessions{exp_nber,'sub'}{:} '_ses-' sessions{exp_nber,'ses'}{:}...
        '_LFP_pvalue_trial_table_refLaplacian.png']);
    if exist(dir_path,"dir")
        save(file_path,'-struct','p_values');
        saveas(gcf,fig_path);
        disp("Results saved as .mat and .png");
        disp("--------------------");
        clear PSD_results p_values pvalue_trial p_values_comp;
        close all
        root_dir='Y:\';
        [~,sessions]=fetch_flicker_subjectIDs(root_dir,'flickerneuro');
        p_values.ses = join([sessions.sub,sessions.ses],'_',2);
    else
        disp("File path not found!")
    end
   

end

% figure
% plot_PSD('80Hz-AV','1Ld5-1Ld4/1Ld6',PSD_results,PSD_results.label,PSD_results.condition,condition_color('80Hz-AV'),0,1,1)% the 0 is for std dev 
% hold on;
% plot_PSD('Baseline','1Ld5-1Ld4/1Ld6',PSD_results,PSD_results.label,PSD_results.condition,[0 0 0],1,0,1)
% hold off;
       