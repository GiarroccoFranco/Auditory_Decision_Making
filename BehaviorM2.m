clearvars;
close all;
cd('C:\Users\giarroccof2\OneDrive - National Institutes of Health\Franco\Corrie Data\Task') % Folder with data



load M2_BehavioralData
load M2_BehavioralCodes
load M2_A1_Data_100ms
load M2_PFC_Data_100ms


for s=1: size (M2_A1_Data_100ms,2)

    clearvars -except p_NonMatch p_Match ...
        Performance_Del M1Perf_MNM M1Perf_del ...
        s p_del Acc_PFC Acc_AC M2_A1_Data_100ms mean_RTs_Match_Corr mean_RTs_Match NtrialsSessions PFC_Data AC_Data  Eig_Values_AC_Sess_Match  Eig_Values_AC_Sess Eig_Values_PFC_Sess Eig_Values_PFC_Sess_NonMatch ...
        nk k Times Smoothvalue SelectedDelays Times_fig M2_A1_Data_100ms  M2_PFC_Data_100ms n_neu_across_sess_Mega  M2_BehavioralCodes M2_BehavioralData T_M  T_NonM  CRTR CLTL CRTR CRTL Monkey_acc_sessions


    size (M2_A1_Data_100ms,2) - s

    An_AC_Sess    = M2_A1_Data_100ms{s};
    An_PFC_Sess   = M2_PFC_Data_100ms{s};

    An_Behav_Sess = M2_BehavioralData{s};

    All_Trials    = 1:size(An_AC_Sess,1);

    catchtr       = M2_BehavioralCodes{s}.iscatch;


    Catch_Trials  = find (catchtr==1);

    acc           = M2_BehavioralCodes{s}.acc;
    acc_Nan       = find(isnan(acc));

    Catch_and_Acc_nan=intersect(Catch_Trials,acc_Nan);
    Acc_Nan_non_Catch=setdiff(acc_Nan,Catch_Trials);

    trials_2_check=randsample(Catch_and_Acc_nan, 10);

    checkTrials_Markers = M2_BehavioralCodes{s}.t2rt  (Acc_Nan_non_Catch );
    checkTrials_behav         = An_Behav_Sess.A(Acc_Nan_non_Catch,:)
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




    Match_Trials_2=zeros(size(TargetSide));
    Match_2         = TargetSide==CueLocation;
    idx_Match_2     = find(Match_2==1);
    idx_nonMatch_2  = find(Match_2==0);
    Match_Trials_2(idx_Match_2)=1;
    Match_Trials_2(idx_nonMatch_2)=-1;


    MnM_Acc_Nan_non_Catch= Match_Trials_2(Acc_Nan_non_Catch);

    Match_Corr_Trials  = zeros(size(TargetSide_Correct_trials));
    Match_Corr         = CueSide_Correct_trials==TargetSide_Correct_trials;
    idx_Match_Corr     = find(Match_Corr==1);
    idx_nonMatch_Corr  = find(Match_Corr==0);
    Match_Corr_Trials    (idx_Match_Corr)= 1;
    Match_Corr_Trials    (idx_nonMatch_Corr)= -1;

    Match_Trials  = zeros(size(TargetSide_Used_trials));
    Match         = CueSide_Used_trials==TargetSide_Used_trials;
    idx_match     = find(Match==1);
    idx_nonmatch  = find(Match==0);
    Match_Trials    (idx_match)= 1;
    Match_Trials    (idx_nonmatch)= -1;

    AllDelayTimes = An_Behav_Sess.A(:,8);

    MnM_Acc_Nan_non_Catch= Match_Trials_2(Acc_Nan_non_Catch);



    Del_1 = find(AllDelayTimes==-1000);
    Del_2 = find(AllDelayTimes==-1300);
    Del_3 = find(AllDelayTimes==-1800);

    Delay_1 = intersect(Del_1,Used_Trials);
    Delay_2 = intersect(Del_2,Used_Trials);
    Delay_3 = intersect(Del_3,Used_Trials);

    Delay_1_2nd = intersect(Del_1,Acc_Nan_non_Catch);
    Delay_2_2nd = intersect(Del_2,Acc_Nan_non_Catch);
    Delay_3_2nd = intersect(Del_3,Acc_Nan_non_Catch);

    Delay_All = union(union(Delay_1,Delay_2),Delay_3);

    Delay_Corr_1 = intersect(Del_1,acc_1);
    Delay_Corr_2 = intersect(Del_2,acc_1);
    Delay_Corr_3 = intersect(Del_3,acc_1);
    Delay_Corr_All = union(union(Delay_Corr_1,Delay_Corr_2),Delay_Corr_3);

    Del_used=AllDelayTimes(Used_Trials);
    Corr_Del_used=AllDelayTimes(acc_1);

    Del_with_Nan = AllDelayTimes(Acc_Nan_non_Catch);



    %%%% for used trials
    Del_Match=[Del_used Match_Trials] ;
    idxDel1=find(Del_Match(:,1)==-1000);
    idxDel2=find(Del_Match(:,1)==-1300);
    idxDel3=find(Del_Match(:,1)==-1800);


    %%%% for Those nan trials
    Del_Match_nan=[Del_with_Nan MnM_Acc_Nan_non_Catch] ;
    idxDel1nan=find(Del_Match_nan(:,1)==-1000);
    idxDel2nan=find(Del_Match_nan(:,1)==-1300);
    idxDel3nan=find(Del_Match_nan(:,1)==-1800);
    %%%% find Match and Non-match for Those nan trials
    idxDel_Mnan=find(Del_Match_nan(:,2)==1);
    idxDel_NMnan=find(Del_Match_nan(:,2)==-1);
    %%%% find Match and Non-match for del 1 2 and 3 for Those nan trials
    del_1_Mnan=intersect(idxDel1nan,idxDel_Mnan);
    del_1_NMnan=intersect(idxDel1nan,idxDel_NMnan);
    del_2_Mnan=intersect(idxDel2nan,idxDel_Mnan);
    del_2_NMnan=intersect(idxDel2nan,idxDel_NMnan);
    del_3_Mnan=intersect(idxDel3nan,idxDel_Mnan);
    del_3_NMnan=intersect(idxDel3nan,idxDel_NMnan);



    %%%% find Match and Non-match
    idxDel_M=find(Del_Match(:,2)==1);
    idxDel_NM=find(Del_Match(:,2)==-1);

    %%%% find Match and Non-match for del 1 2 and 3
    del_1_M=intersect(idxDel1,idxDel_M);
    del_1_NM=intersect(idxDel1,idxDel_NM);
    del_2_M=intersect(idxDel2,idxDel_M);
    del_2_NM=intersect(idxDel2,idxDel_NM);
    del_3_M=intersect(idxDel3,idxDel_M);
    del_3_NM=intersect(idxDel3,idxDel_NM);


    %%%% for correct trials
    Corr_Del_Match=[Corr_Del_used Match_Corr_Trials] ;
    Corr_idxDel1=find(Corr_Del_Match(:,1)==-1000);
    Corr_idxDel2=find(Corr_Del_Match(:,1)==-1300);
    Corr_idxDel3=find(Corr_Del_Match(:,1)==-1800);


    %%%% find Match and Non-match
    Corr_idxDel_M=find(Corr_Del_Match(:,2)==1);
    Corr_idxDel_NM=find(Corr_Del_Match(:,2)==-1);
    %%%% find Correct Match and Non-match for del 1 2 and 3
    Corr_del_1_M=intersect(Corr_idxDel1,Corr_idxDel_M);
    Corr_del_1_NM=intersect(Corr_idxDel1,Corr_idxDel_NM);
    Corr_del_2_M=intersect(Corr_idxDel2,Corr_idxDel_M);
    Corr_del_2_NM=intersect(Corr_idxDel2,Corr_idxDel_NM);
    Corr_del_3_M=intersect(Corr_idxDel3,Corr_idxDel_M);
    Corr_del_3_NM=intersect(Corr_idxDel3,Corr_idxDel_NM);

    % p_Match(s,:)= [ numel(Corr_del_1_M) numel(Corr_del_2_M)  numel(Corr_del_3_M)]./ [numel(del_1_M)+numel(del_1_Mnan) numel(del_2_M)+numel(del_2_Mnan) numel(del_3_M)+numel(del_3_Mnan)];
    % p_NonMatch(s,:)= [ numel(Corr_del_1_NM) numel(Corr_del_2_NM)  numel(Corr_del_3_NM)]./ [numel(del_1_NM)+numel(del_1_NMnan) numel(del_2_NM)+numel(del_2_NMnan) numel(del_3_NM)+numel(del_3_NMnan)];

    p_Match(s,:)= [ numel(Corr_del_1_M) numel(Corr_del_2_M)  numel(Corr_del_3_M)]./ [numel(del_1_M) numel(del_2_M) numel(del_3_M)];
    p_NonMatch(s,:)= [ numel(Corr_del_1_NM) numel(Corr_del_2_NM)  numel(Corr_del_3_NM)]./ [numel(del_1_NM) numel(del_2_NM) numel(del_3_NM)];



    p_del(s,:)=[numel(Delay_Corr_1) numel(Delay_Corr_2) numel(Delay_Corr_3)] ./[numel(Delay_1)+numel(Delay_1_2nd) numel(Delay_2)+numel(Delay_2_2nd) numel(Delay_3)+numel(Delay_3_2nd)];











    RTs_1= M2_BehavioralCodes{s}.t1rt  ;
    Rts_used_Trials=RTs_1(Used_Trials);
    Rts_Match_Trials = Rts_used_Trials(idxDel_M);
    mean_RTs_Match(s,:)=nanmean(Rts_Match_Trials);

    Rts_corr_Trials=RTs_1(acc_1);
    Rts_Match_corrTrials = Rts_corr_Trials(Corr_idxDel_M);
    mean_RTs_Match_Corr(s,:)=nanmean(Rts_Match_corrTrials);






end
PerformaceMonkey=[p_Match p_NonMatch];

% for RTs

figure('Units', 'pixels', 'Position', [500 500 350 250]);
bar(mean(mean_RTs_Match));hold on
errorbar(mean(mean_RTs_Match),std(mean_RTs_Match)/sqrt(numel(mean_RTs_Match)));
xlim([0 7])
ylim([0 .7])
ylabel ('Response time (s)')




% for performance Match vs Non-Match
data=PerformaceMonkey;
% Suppose data is n x 6
n = size(data,1);
meanData = mean(data, 1);                  % 1 x 6
semData  = std(data, 0, 1) ;      % 1 x 6
% Arrange means and SEMs into 2 x 3
M = [    meanData(1), meanData(4), meanData(2)   meanData(5), meanData(3), meanData(6)];
E = [   semData(1), semData(4), semData(2) semData(5), semData(3), semData(6)];

Mt = M;  % 3 x 2
Et = E';  % 3 x 2
x_bar= [1 2 3.5 4.5 6 7 ];
figure('Units', 'pixels', 'Position', [500 500 350 250]);
b = bar(x_bar,Mt,'grouped');
hold on;
% Add error bars
errorbar(x_bar, M, E, ...
    'k', 'linestyle', 'none', 'LineWidth', 1);
ylabel ('Accuracy')

hold off;