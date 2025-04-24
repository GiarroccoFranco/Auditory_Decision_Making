clearvars
cd('C:\Users\giarroccof2\OneDrive - National Institutes of Health\Franco\Corrie Data\Task')


load M1_PFC_Data_50ms_Cue.mat
load M1_A1_Data_50ms_Cue.mat
load M1_BehavioralData.mat
load M1_BehavioralCodes.mat

Sensory_Window=20:31;
Times=-1000:50:3500;
Times_fig=-900:50:3500;




for s=1: size (M1_A1_Data_50ms_Cue,2)
    clearvars -except M1_PFC_Data_50ms_Cue M1_A1_Data_50ms_Cue Sensory_Window ...
        Times Times_fig s  R_s_Predictions M1_BehavioralData M1_BehavioralCodes

    size (M1_A1_Data_50ms_Cue,2) - s

    An_AC_Sess    = M1_A1_Data_50ms_Cue{s};
    An_PFC_Sess   = M1_PFC_Data_50ms_Cue{s};

    An_Behav_Sess = M1_BehavioralData{s};

    All_Trials    = 1:size(An_AC_Sess,1);

    catchtr       = M1_BehavioralCodes{s}.iscatch;
    Catch_Trials  = find (catchtr==1);

    acc           = M1_BehavioralCodes{s}.acc;
    acc_Nan       = find(isnan(acc));
    RemovedTrials=[];
    RemovedTrials = union(Catch_Trials,acc_Nan);
    All_Trials(RemovedTrials)   = [];

    Used_Trials   = All_Trials;

    pos_acc=find(acc==1);
    acc_1= intersect(Used_Trials,pos_acc);

    neg_acc=find(acc==0);
    acc_0= intersect(Used_Trials,neg_acc);

    RTs_1= M1_BehavioralCodes{s}.t1rt  ;

    RTs_1_Corr =RTs_1(Used_Trials);

    p33 = prctile(RTs_1_Corr(:,1), 33);
    p66 = prctile(RTs_1_Corr(:,1), 66);

    Trials_RT = nan(size(RTs_1_Corr,1),3);


    % Categorize the trials based on RT
    for i = 1:size(RTs_1_Corr,1)

        if RTs_1_Corr(i,1) <= p33
            Trials_RT(i,1) =  RTs_1_Corr(i,1); % Collect trial numbers for short RTs
        elseif RTs_1_Corr(i,1) > p33 && RTs_1_Corr(i,1) <= p66
            Trials_RT(i,2) =  RTs_1_Corr(i,1); % Collect trial numbers for short RTs
        else
            Trials_RT(i,3) =  RTs_1_Corr(i,1); % Collect trial numbers for short RTs
        end
    end

    acc_RTs1=Used_Trials(find(isnan(Trials_RT(:,1))==0));
    acc_RTs2=Used_Trials(find(isnan(Trials_RT(:,2))==0));
    acc_RTs3=Used_Trials(find(isnan(Trials_RT(:,3))==0));

    acc_RTs={acc_RTs1,acc_RTs2,acc_RTs3};

    CueLocation                 = An_Behav_Sess.A(:,7);
    CueLocation(CueLocation==2) = -1;
    CueSide_Used_trials         = CueLocation(Used_Trials);
    CueSide_Correct_trials      = CueLocation(acc_1);

    TargetSide    = M1_BehavioralCodes{s}.tside;
    TargetSide(TargetSide==2) = -1;
    TargetSide_Used_trials         = TargetSide(Used_Trials);
    TargetSide_Correct_trials      = TargetSide(acc_1);

    Match_Trials  = zeros(size(TargetSide_Used_trials));
    Match         = CueSide_Used_trials==TargetSide_Used_trials;
    idx_match     = find(Match==1);
    idx_nonmatch  = find(Match==0);
    Match_Trials    (idx_match)= 1;
    Match_Trials    (idx_nonmatch)= -1;

    Match_Corr_Trials  = zeros(size(TargetSide_Correct_trials));
    Match_Corr         = CueSide_Correct_trials==TargetSide_Correct_trials;
    idx_Match_Corr     = find(Match_Corr==1);
    idx_nonMatch_Corr  = find(Match_Corr==0);
    Match_Corr_Trials    (idx_Match_Corr)= 1;
    Match_Corr_Trials    (idx_nonMatch_Corr)= -1;

    AllDelayTimes = An_Behav_Sess.A(Used_Trials,8);

    Delay_1=  find(AllDelayTimes==-1000);
    Delay_2=  find(AllDelayTimes==-1300);
    Delay_3=  find(AllDelayTimes==-1800);

    AllDelay_CorrTimes = An_Behav_Sess.A(acc_1,8);

    Delay_Corr_1=  find(AllDelay_CorrTimes==-1000);
    Delay_Corr_2=  find(AllDelay_CorrTimes==-1300);
    Delay_Corr_3=  find(AllDelay_CorrTimes==-1800);

    RTs_1_Used_Trials =RTs_1(Used_Trials);

    RT_not_nan=find(isnan(RTs_1_Used_Trials)==0);
    Used_Trials_with_RTs=Used_Trials(RT_not_nan);
    RT_values_UsedTrials =  RTs_1_Used_Trials(RT_not_nan);

    median_Used_Trials = median( RTs_1_Used_Trials, "omitnan");
    Slow_RTidx=find (RTs_1>median_Used_Trials);
    Fast_RTidx=find (RTs_1<median_Used_Trials);

    SlowRT_Used_Trials= intersect(Slow_RTidx,Used_Trials);
    FastRT_Used_Trials= intersect(Fast_RTidx,Used_Trials);
    Slow_RT=RTs_1 (SlowRT_Used_Trials);
    MnM_Label_SlowRT=Match_Trials(find(ismember(Used_Trials,SlowRT_Used_Trials)));
    MnM_Label_FastRT=Match_Trials(find(ismember(Used_Trials,FastRT_Used_Trials)));

    trials_For_Regression=Used_Trials;

    BetaCue_AC=[]; BetaTg_AC=[]; NormCue_AC=[]; NormTg_AC=[];
%% First regression used to define 1D-dimension
    % A1 regression
    for AC_neu= 1 : size (An_AC_Sess,3) % # neurons
        AC_single_Neu=[]; AC_single_Neu = An_AC_Sess(:,:,AC_neu);

        for bin=1:size(AC_single_Neu,2)
            B = fitlm([CueSide_Used_trials],AC_single_Neu(trials_For_Regression,bin));

            BetaCue_AC(AC_neu,bin)=B.Coefficients{2,1};
        end
    end

    BetaCue_PFC=[]; BetaTg_PFC=[]; NormCue_PFC=[]; NormTg_PFC=[];
    % PFC regression
    for PFC_neu= 1 : size (An_PFC_Sess,3) % # neurons
        PFC_single_Neu=[]; PFC_single_Neu = An_PFC_Sess(:,:,PFC_neu);

        for bin=1:size(PFC_single_Neu,2)
            B = fitlm([CueSide_Used_trials],PFC_single_Neu(trials_For_Regression,bin));

            BetaCue_PFC(PFC_neu,bin)=B.Coefficients{2,1};
        end
    end


    %% Define 1D-dimension
    NormCue_AC=vecnorm(BetaCue_AC);
    Max_Norm_Cue(s,1)=max(NormCue_AC(Sensory_Window));
    bin_betaCue_AC=find(NormCue_AC==max(NormCue_AC(Sensory_Window)));
    BetaCue_Vector_AC=BetaCue_AC(:,bin_betaCue_AC(1));

    NormCue_PFC=vecnorm(BetaCue_PFC);
    bin_betaCue_PFC=find(NormCue_PFC==max(NormCue_PFC(Sensory_Window)));
    BetaCue_Vector_PFC=BetaCue_PFC(:,bin_betaCue_PFC(1));  %?
    Max_Norm_Cue(s,2)=max(NormCue_PFC(Sensory_Window));

    betas(s,1)=Times(bin_betaCue_AC);
    betas(s,2)=Times(bin_betaCue_PFC);

    AC_PFC_Data= {permute(An_AC_Sess, [3 2 1]), permute(An_PFC_Sess, [3 2 1])} ;
    IT_Micro_Tg_MnM_Values_All_Trials_2= [Max_Norm_Cue betas ];
    BetaCue_Vectors= {BetaCue_Vector_AC,BetaCue_Vector_PFC };
    % NormTg_AC=vecnorm(BetaTg_AC);
    % bin_betaTg_AC=find(NormTg_AC==max(NormTg_AC(Sensory_Window)));
    % BetaTg_Vector_AC=BetaTg_AC(:,bin_betaTg_AC(1));


    %% Targeted dimensionality reduction 
    n_areas=2;
    
    Proj_B_Cue=nan( size(trials_For_Regression,2), size(AC_PFC_Data{1,1},2),    size(AC_PFC_Data,2)); % Trials x bin x Areas
    for n_Areas= 1:n_areas
        PermData=[]; PermData=AC_PFC_Data{1,n_Areas};
        for     Trial_Beta=1:size(trials_For_Regression,2)
            for tb=1:size (PermData,2)
                Proj_B_Cue(Trial_Beta,tb,n_Areas)=dot(PermData(:,tb,trials_For_Regression(Trial_Beta)),BetaCue_Vectors{1,n_Areas});
            end
        end
    end


    %% To remove same area intrinsic temporal dependencies within each area 
    Residual_SameArea=[];
    for  n_Areas= 1:n_areas
        Area_Feature= Proj_B_Cue(:,:,n_Areas);

        for b= 1:size(Proj_B_Cue,2)-1
            Self_Prediction=fitlm(Area_Feature(:,b),Area_Feature(:,b+1));
            Residual_SameArea(:,b,n_Areas)=Self_Prediction.Residuals{:,2}  ;
        end
    end

    % %%%%%
    % %%%%% to validate the previous autocorrelation model...
    % %%%%% 
    % %%%%%

    % residuals = Residual_SameArea(:,:,1);
    %
    % max_lag = 5; %%%%% Define the maximum lag to test
    % num_trials = size(residuals, 1); %%%%% Number of trials
    % num_lags = 2 * max_lag + 1; %%%%%% Lags from -max_lag to max_lag
    % xcorr_vals = zeros(num_trials, num_lags); %%%%%% Store cross-correlations
    %
    % for trial = 1:num_trials
    %     xcorr_vals(trial, :) = xcorr(residuals(trial, :), max_lag, 'coeff'); %%%%%%% Compute xcorr for each trial
    % end
    %
    % %%%%%%%% Average across trials to obtain session-level cross-correlation
    % mean_xcorr = mean(xcorr_vals, 1);
    %
    % %%%%%%%% Plot results
    % lags = -max_lag:max_lag;
    %
    %
    % %%%%%%%%% this is to plot the autocorrelation to ensure that the previous
    % %%%%%%%%% autoregression whitened the time series
    % figure;
    % stem(lags,mean_xcorr)
    % xlabel('Lag');
    % ylabel('Autocorrelation coefficient');
    % hold off;
    %
    % %%%%%%% residuals: matrix of size (nTrials x nTimeBins)
    % %%%%%%% Each row is one trial's residual time series
    %
    % maxLag = 5;                  %%%%%%%% Maximum lag to evaluate
    % numTrials = size(residuals, 1);
    %
    % %%%%%%%% Preallocate to store the xcorr results for each trial
    % ccVals = zeros(numTrials, 2*maxLag + 1);
    %
    % for tr = 1:numTrials
    %     %%%%%%%% Compute xcorr for one trial, up to maxLag, normalized by 'coeff'
    %     [cc, lags] = xcorr(residuals(tr, :), maxLag, 'coeff');
    %     ccVals(tr, :) = cc;     %%%%%%%% Store the cross-correlation values
    % end
    %
    % %%%%%%% Average across trials
    % meanCC = mean(ccVals, 1);
    % semCC  = std(ccVals, [], 1) / sqrt(numTrials);  %%%%%%%% Optional SEM
    %
    % %%%%%%%%%% Plot the mean cross-correlation
    % figure;
    % errorbar(lags, meanCC, semCC, '-o', 'LineWidth', 2, 'MarkerFaceColor','black');
    % xlabel('Lag');
    % ylabel('Autocorrelation coefficient');
    % grid on;
    %
    % alpha = 0.05;
    % for i = 1:length(lags)
    %     if lags(i) ~= 0
    %         [~, p] = ttest(ccVals(:, i), 0, 'Alpha', alpha);
    %         fprintf('Lag %d -> p = %.4f\n', lags(i), p);
    %     end
    % end

    %%%%%%%%%
    %%%%%%%%% end of autocorrelation validation
    %%%%%%%%%

    %% Perform Cross-area prediction
    for  n_Areas= 1:n_areas

        if n_Areas==1
            TargetAreas=1; SourceAreas=2;
        elseif n_Areas==2
            TargetAreas=2; SourceAreas=1;
        end

        for b= 1:size(Residual_SameArea,2)-1
            IT_Prediction=fitlm(Residual_SameArea(:,b,SourceAreas),Residual_SameArea(:,b+1,TargetAreas));
            R_s_Predictions(s,b,n_Areas)=IT_Prediction.Rsquared.Adjusted;
        end
    end

end
options.x_axis = Times_fig(1:end) ;
Colors =  [ [250 160 10]/255;  [4 170 230]/255; ];
options.s      = 2;
options.color_PFC = Colors(1,:);
options.color_AC  = Colors(2,:);



figure('Units', 'pixels', 'Position', [500 500 350 200 ]);
plot_areaerrorbarAC(  R_s_Predictions(:,:,2),options); hold on
plot_areaerrorbarPFC( R_s_Predictions(:,:,1),options); hold on
% plot_areaerrorbarAC(  R_s_Predictions([1:15 17:end],:,2),options); hold on
% plot_areaerrorbarPFC( R_s_Predictions([1:15 17:end],:,1),options); hold on
line([0 0], [-.01 0.04],'LineStyle','--','color','k')
xlim([-500 1000])
ylim([-.01 .04])
xlabel('Time from Target (ms)')
ylabel('r^2')

