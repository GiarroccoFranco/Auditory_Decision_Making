clearvars;
% close all;
cd('C:\Users\giarroccof2\OneDrive - National Institutes of Health\Franco\Corrie Data\Task')


load n_neurons_M2
load M2_BehavioralData
load M2_BehavioralCodes
load M2_A1_Data_100ms
load M2_PFC_Data_100ms



SelectedDelays= 3 ;
Times=-1000:25:3500;% for cue alignment 20 ms

AC_Data = nan(sum(n_neurons_M2(1,:)),size(M2_PFC_Data_100ms{1, 1}  ,2)*12) ;  PFC_Data = nan(sum(n_neurons_M2(2,:)),size(M2_PFC_Data_100ms{1, 1}  ,2)*12) ;


for s=1: size (M2_A1_Data_100ms,2)

    clearvars -except s Acc_PFC Acc_AC M2_A1_Data_100ms NtrialsSessions PFC_Data AC_Data  Eig_Values_AC_Sess_Match  Eig_Values_AC_Sess Eig_Values_PFC_Sess Eig_Values_PFC_Sess_NonMatch ...
        nk k Times Smoothvalue SelectedDelays Times_fig M2_A1_Data_100ms  M2_PFC_Data_100ms n_neurons_M2  M2_BehavioralCodes M2_BehavioralData T_M  T_NonM  CRTR CLTL CRTR CRTL Monkey_acc_sessions


    size (M2_A1_Data_100ms,2) - s

    An_AC_Sess    = M2_A1_Data_100ms{s};
    An_PFC_Sess   = M2_PFC_Data_100ms{s};

    An_Behav_Sess = M2_BehavioralData{s};

    All_Trials    = 1:size(An_AC_Sess,1);

    catchtr       = M2_BehavioralCodes{s}.iscatch;
    Catch_Trials  = find (catchtr==1);

    acc           = M2_BehavioralCodes{s}.acc;
    acc_Nan       = find(isnan(acc));

    RemovedTrials = union(Catch_Trials,acc_Nan);
    All_Trials(RemovedTrials)   = [];

    Used_Trials   = All_Trials;

    pos_acc=find(acc==1);
    acc_1= intersect(Used_Trials,pos_acc);

    neg_acc=find(acc==0);
    acc_0= intersect(Used_Trials,neg_acc);


    CueLocation                 = An_Behav_Sess.A(:,7);
    CueLocation(CueLocation==2) = -1;
    CueSide_Used_trials         = CueLocation(Used_Trials);
    CueSide_Correct_trials      = CueLocation(acc_1);


    TargetSide    = M2_BehavioralCodes{s}.tside;
    TargetSide(TargetSide==2) = -1;
    TargetSide_Used_trials         = TargetSide(Used_Trials);
    TargetSide_Correct_trials      = TargetSide(acc_1);


    Match_Trials  = zeros(size(TargetSide_Used_trials));
    Match         = CueSide_Used_trials==TargetSide_Used_trials;
    idx_match     = find(Match==1);
    idx_nonmatch  = find(Match==0);
    Match_Trials    (idx_match)= 1;
    Match_Trials    (idx_nonmatch)= -1;


    Match_Trials  = zeros(size(An_AC_Sess,1),1);
    Match         = CueLocation==TargetSide;
    idx_match     = find(Match==1);
    idx_nonmatch  = find(Match==0);
    Match_Trials    (idx_match)= 1;
    Match_Trials    (idx_nonmatch)= -1;

    AllDelayTimes = An_Behav_Sess.A(:,8);

    Del_1 = find(AllDelayTimes==-1000);
    Del_2 = find(AllDelayTimes==-1300);
    Del_3 = find(AllDelayTimes==-1800);
    Delay_1 = intersect(Del_1,Used_Trials);
    Delay_2 = intersect(Del_2,Used_Trials);
    Delay_3 = intersect(Del_3,Used_Trials);
    Delay_All = union(union(Delay_1,Delay_2),Delay_3);

    Delay_Corr_1 = intersect(Del_1,acc_1);
    Delay_Corr_2 = intersect(Del_2,acc_1);
    Delay_Corr_3 = intersect(Del_3,acc_1);
    Delay_Corr_All = union(union(Delay_Corr_1,Delay_Corr_2),Delay_Corr_3);





    for Del = 1:SelectedDelays
        if SelectedDelays==1
            Delay_s = 0 ;
        elseif SelectedDelays==3
            Delay_s = Del ;
        end

        if Delay_s ==1


            Cue_Del_1=CueLocation(Delay_1);
            Cue_Del_1_Corr=CueLocation(Delay_Corr_1);

            TargetSide_Del_1= TargetSide(Delay_1);
            TargetSide_Del_1_Corr= TargetSide(Delay_Corr_1);

            Match_Del_1=Match_Trials(Delay_1);
            Match_Del_1_Corr=Match_Trials(Delay_Corr_1);

            Cue_Left_Del_1=Delay_Corr_1( find(Cue_Del_1_Corr== -1));
            Cue_Right_Del_1=Delay_Corr_1( find(Cue_Del_1_Corr== 1));

            Target_Left_Del_1=Delay_Corr_1( find(TargetSide_Del_1_Corr== -1));
            Target_Right_Del_1=Delay_Corr_1( find(TargetSide_Del_1_Corr== 1));

            Cue_Left_Target_Match = intersect(Cue_Left_Del_1,Target_Left_Del_1);
            Cue_Left_Target_NonMatch = intersect(Cue_Left_Del_1,Target_Right_Del_1);
            Cue_Right_Target_Match = intersect(Cue_Right_Del_1,Target_Right_Del_1);
            Cue_Right_Target_NonMatch = intersect(Cue_Right_Del_1,Target_Left_Del_1);



        elseif Delay_s ==2
            Cue_Del_2=CueLocation(Delay_2);
            Cue_Del_2_Corr=CueLocation(Delay_Corr_2);

            TargetSide_Del_2= TargetSide(Delay_2);
            TargetSide_Del_2_Corr= TargetSide(Delay_Corr_2);

            Match_Del_2=Match_Trials(Delay_2);
            Match_Del_2_Corr=Match_Trials(Delay_Corr_2);

            Cue_Left_Del_2=Delay_Corr_2( find(Cue_Del_2_Corr== -1));
            Cue_Right_Del_2=Delay_Corr_2( find(Cue_Del_2_Corr== 1));

            Target_Left_Del_2=Delay_Corr_2( find(TargetSide_Del_2_Corr== -1));
            Target_Right_Del_2=Delay_Corr_2( find(TargetSide_Del_2_Corr== 1));

            Cue_Left_Target_Match = intersect(Cue_Left_Del_2,Target_Left_Del_2);
            Cue_Left_Target_NonMatch = intersect(Cue_Left_Del_2,Target_Right_Del_2);
            Cue_Right_Target_Match = intersect(Cue_Right_Del_2,Target_Right_Del_2);
            Cue_Right_Target_NonMatch = intersect(Cue_Right_Del_2,Target_Left_Del_2);


        elseif Delay_s ==3
            Cue_Del_3=CueLocation(Delay_3);
            Cue_Del_3_Corr=CueLocation(Delay_Corr_3);

            TargetSide_Del_3= TargetSide(Delay_3);
            TargetSide_Del_3_Corr= TargetSide(Delay_Corr_3);

            Match_Del_3=Match_Trials(Delay_3);
            Match_Del_3_Corr=Match_Trials(Delay_Corr_3);

            Cue_Left_Del_3=Delay_Corr_3( find(Cue_Del_3_Corr== -1));
            Cue_Right_Del_3=Delay_Corr_3( find(Cue_Del_3_Corr== 1));

            Target_Left_Del_3=Delay_Corr_3( find(TargetSide_Del_3_Corr== -1));
            Target_Right_Del_3=Delay_Corr_3( find(TargetSide_Del_3_Corr== 1));

            Cue_Left_Target_Match = intersect(Cue_Left_Del_3,Target_Left_Del_3);
            Cue_Left_Target_NonMatch = intersect(Cue_Left_Del_3,Target_Right_Del_3);
            Cue_Right_Target_Match = intersect(Cue_Right_Del_3,Target_Right_Del_3);
            Cue_Right_Target_NonMatch = intersect(Cue_Right_Del_3,Target_Left_Del_3);



        end





        Trials_ToAnalize={};
        Trials_ToAnalize={Cue_Left_Target_Match,Cue_Left_Target_NonMatch,Cue_Right_Target_Match,Cue_Right_Target_NonMatch};





        AC_C_L_T_L=[];   AC_C_L_T_R=[];   AC_C_R_T_R=[];  AC_C_R_T_L =[];    AC_Wr_Trials_L =[];

        % AC regression
        for AC_neu= 1 : size (An_AC_Sess,3) % # neurons

            AC_single_Neu=[]; AC_single_Neu = An_AC_Sess(:,:,AC_neu);


            AC_C_L_T_L (AC_neu,:) = mean(AC_single_Neu(Trials_ToAnalize{1,1},:));
            AC_C_L_T_R (AC_neu,:) = mean(AC_single_Neu(Trials_ToAnalize{1,2},:));

            AC_C_R_T_R (AC_neu,:) = mean(AC_single_Neu(Trials_ToAnalize{1,3},:));
            AC_C_R_T_L (AC_neu,:) = mean(AC_single_Neu(Trials_ToAnalize{1,4},:));


        end


        M_AC(:,:,Del)=horzcat (AC_C_L_T_L, AC_C_L_T_R, AC_C_R_T_R, AC_C_R_T_L);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % # of components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        PFC_C_L_T_L=[];   PFC_C_L_T_R=[];   PFC_C_R_T_R=[];  PFC_C_R_T_L =[];    PFC_Wr_Trials_L =[];

        for PFC_neu= 1 : size (An_PFC_Sess,3) % # neurons

            PFC_single_Neu=[]; PFC_single_Neu = An_PFC_Sess(:,:,PFC_neu);


            PFC_C_L_T_L (PFC_neu,:) = mean(PFC_single_Neu(Trials_ToAnalize{1,1},:));
            PFC_C_L_T_R (PFC_neu,:) = mean(PFC_single_Neu(Trials_ToAnalize{1,2},:));

            PFC_C_R_T_R (PFC_neu,:) = mean(PFC_single_Neu(Trials_ToAnalize{1,3},:));
            PFC_C_R_T_L (PFC_neu,:) = mean(PFC_single_Neu(Trials_ToAnalize{1,4},:));



        end


        M_PFC(:,:,Del)=horzcat (PFC_C_L_T_L, PFC_C_L_T_R, PFC_C_R_T_R, PFC_C_R_T_L);



    end



    %%
  

    M_AC_all_Del = [];    M_PFC_all_Del = [];

    if Delay_s~=0
        M_AC_all_Del=horzcat(M_AC(:,:,1) , M_AC(:,:,2), M_AC(:,:,3) ) ;
        M_PFC_all_Del=horzcat(M_PFC(:,:,1) , M_PFC(:,:,2), M_PFC(:,:,3) ) ;
    else
        M_AC_all_Del =M_AC;
        M_PFC_all_Del=M_PFC;
    end
    M_AC=[]; M_PFC=[];
    if s ==1 && Delay_s~=0
        AC_Data(1:n_neurons_M2(1,1),:)=M_AC_all_Del;
        PFC_Data(1:n_neurons_M2(2,1),:)=M_PFC_all_Del;
    elseif s >1 && Delay_s~=0
        AC_Data(sum(n_neurons_M2(1,1:s-1))+1:sum(n_neurons_M2(1,1:s)),:)=M_AC_all_Del;
        PFC_Data(sum(n_neurons_M2(2,1:s-1))+1:sum(n_neurons_M2(2,1:s)),:)=M_PFC_all_Del;
    end
    if s ==1 && Delay_s==0
        AC_Data(1:n_neurons_M2(1,1),:)=M_AC_all_Del;
        PFC_Data(1:n_neurons_M2(2,1),:)=M_PFC_all_Del;
    elseif s >1 && Delay_s==0
        AC_Data(sum(n_neurons_M2(1,1:s-1))+1:sum(n_neurons_M2(1,1:s)),:)=M_AC_all_Del;
        PFC_Data(sum(n_neurons_M2(2,1:s-1))+1:sum(n_neurons_M2(2,1:s)),:)=M_PFC_all_Del;
    end


end




%%

% close all
k = 3; % Number of components to retain
% AC
AC_Sessions=[];score_AC=[];
AC_Sessions=AC_Data;


[coeff,score_ac,latent,tsquared,explained_ac,mu] = pca(AC_Sessions');

score_AC=smoothdata(score_ac(:,1:k ),1,"gaussian",5);

% - PFC
PFC_Sessions=[]; score_PFC=[];
PFC_Sessions=PFC_Data;

[coeff,score_pfc,latent,tsquared,explained_pfc,mu] = pca(PFC_Sessions');


score_PFC=smoothdata(score_pfc(:,1:k ),1,"gaussian",5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Trajectories 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numConditions    = 12;
binsPerCondition = length(Times);
numFigures  = 3;
condsPerFig = 4;




bin_duration = 25;
timeStart = 10;
cueBin = 41;
targetBins  = [81, 93, 113];

timeEnds = [110, 120, 140];

Ax_end=[1550 1850 2400] ;

delayLabels = ["Short Delay", "Medium Delay", "Long Delay"];

Colors3D = [
    250 160  10
    250  70  10
    4 170 230
    10  90 170
    ] / 255;

greyColor = [0.5 0.5 0.5];


conditionLabels = ["R-R", "R-L", "L-L", "L-R", "5", "6"];

for Areas = 1:2

    if Areas == 1
        data_AC_PFC = score_AC;
        areaName    = 'AC';
    else
        data_AC_PFC = score_PFC;
        areaName    = 'PFC';
    end

    for figIdx = 1:numFigures


delayLabels = ["Short Delay", "Medium Delay", "Long Delay"];

        timeEnd  = timeEnds(figIdx);      % one of [140,160,190]
        timeRange = timeStart:timeEnd;    %
        n_bins    = length(timeRange);    %


        cueIdx = cueBin - timeStart + 1;
        xtimes = ((1:n_bins) - cueIdx) * bin_duration;

        % Conditions
        condStart    = (figIdx-1)*condsPerFig + 1;
        condEnd      = figIdx*condsPerFig;
        currentConds = condStart : condEnd;
        num_trials   = length(currentConds);

        thisTargetBin = targetBins(figIdx);

       
        %======================================================================
        % Euclidean Distance
        %======================================================================
     figure('Units', 'pixels', 'Position', [500 100 350 250],'Name',...
            sprintf('%s: %s', areaName, delayLabels(figIdx)),  'NumberTitle','off');
        hold on;

        data_AC_PFC_forFig = zeros(n_bins * num_trials, 3);

        for iCon = 1:num_trials
            cond     = currentConds(iCon);
            idxStart = (cond-1)*binsPerCondition + 1;
            rowsToUse = idxStart + (timeRange - 1);

            rowOutStart = (iCon-1)*n_bins + 1;
            rowOutEnd   = iCon*n_bins;
            data_AC_PFC_forFig(rowOutStart:rowOutEnd, :) = data_AC_PFC(rowsToUse, :);
        end




        baseline_duration = 500;            % ms   300 for PFC
        baseline_bins     = baseline_duration / bin_duration;  % 25 bins
        baseline_time     = 1:baseline_bins; % first 25 bins as baseline

        pairwise_ED_time = zeros(n_bins, num_trials, num_trials);

        %  Euclidean distances
        for trial1 = 1:num_trials
            for trial2 = trial1:num_trials
                s1 = (trial1 - 1)*n_bins + 1;
                e1 = trial1*n_bins;
                traj1 = data_AC_PFC_forFig(s1:e1, :);

                s2 = (trial2 - 1)*n_bins + 1;
                e2 = trial2*n_bins;
                traj2 = data_AC_PFC_forFig(s2:e2, :);

                distances = sqrt(sum((traj1 - traj2).^2, 2));
                pairwise_ED_time(:, trial1, trial2) = distances;
                pairwise_ED_time(:, trial2, trial1) = distances;
            end
        end


        smoothed_ED_time = smoothdata(pairwise_ED_time, 1, 'gaussian', 1);

        all_baseline_ED = [];
        for trial1 = 1:num_trials
            for trial2 = trial1:num_trials
                baseline_data = smoothed_ED_time(baseline_time, trial1, trial2);
                all_baseline_ED = [all_baseline_ED; baseline_data(:)];
            end
        end

        baseline_ED  = mean(all_baseline_ED);
        baseline_std = std(all_baseline_ED);

        threshold = baseline_ED + 3.5 *baseline_std;

        significant_times = smoothed_ED_time > threshold;

        trial_labels = conditionLabels;

        ColorsED = [
            [250, 160,  10]/255
            [250,  70,  10]/255
            [  4, 170, 230]/255
            [ 10,  90, 170]/255
            [ 50,  50,  50]/255
            [200,  50,  50]/255
            [  0, 200, 150]/255
            [150, 150, 200]/255
            [100, 100,  50]/255
            [250, 200,  50]/255
            [  0,   0,   0]/255
            [  0,   0, 200]/255
            [200,   0,   0]/255
            [  0, 200,   0]/255
            [ 30, 100, 150]/255
            ];

        legend_entries = {};
        pair_idx       = 1;

        for trial1 = 1:num_trials
            for trial2 = trial1+1:num_trials
                legend_entries{pair_idx} = trial_labels(trial1) + "-" + trial_labels(trial2);

                h = plot(xtimes, smoothed_ED_time(:, trial1, trial2), ...
                    'Color', ColorsED(pair_idx, :), 'LineWidth', 1.5);


                set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'on');


                sigMask = significant_times(:, trial1, trial2);
                sigBins = find(sigMask);

                if ~isempty(sigBins)
                    diff_bins = diff(sigBins);
                    breaks    = find(diff_bins > 1);
                    segStart  = [sigBins(1); sigBins(breaks + 1)];
                    segEnd    = [sigBins(breaks); sigBins(end)];

                    for s = 1:length(segStart)
                        thick_line = plot( ...
                            xtimes(segStart(s):segEnd(s)), ...
                            smoothed_ED_time(segStart(s):segEnd(s), trial1, trial2), ...
                            'Color', ColorsED(pair_idx, :), 'LineWidth', 4);

                        set(get(get(thick_line, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
                    end
                end
                pair_idx = pair_idx + 1;
            end
        end

        line([0 0], ylim, 'LineStyle', '--', 'Color', 'k');

        legend(legend_entries, 'Location','best', 'Box','off');

        xlabel('Time from Cue (ms)');
        ylabel('Euclidean Distance ');
        title(sprintf('%s %', areaName, delayLabels(figIdx)));

        axis tight;
        xlim([-500 Ax_end(figIdx)]);
        set(gca, 'FontSize', 12);

        hold off;
    end

end
