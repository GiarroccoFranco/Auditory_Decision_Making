clearvars
cd('C:\Users\giarroccof2\OneDrive - National Institutes of Health\Franco\Corrie Data\Task')


load M1_BehavioralData
load M1_BehavioralCodes
load M1_A1_Data_100ms_Target
load M1_PFC_Data_100ms_Target



SelectedDelays= 4 ;

% Times=-1000:25:3500;% for cue alignment
Times=-2100:25:1700;% for Tg alignment

Acc_PFC=nan(size (M1_A1_Data_100ms_Target,2),size (M1_A1_Data_100ms_Target{1},2),4);
Acc_AC=nan(size (M1_A1_Data_100ms_Target,2),size (M1_A1_Data_100ms_Target{1},2)),4;

% for s= 11
parfor s=1: size (M1_A1_Data_100ms_Target,2)
    % for s= f

    size (M1_A1_Data_100ms_Target,2) - s

    An_AC_Sess    = M1_A1_Data_100ms_Target{s};
    An_PFC_Sess   = M1_PFC_Data_100ms_Target{s};

    An_Behav_Sess = M1_BehavioralData{s};

    All_Trials    = 1:size(An_AC_Sess,1);

    catchtr       = M1_BehavioralCodes{s}.iscatch;
    Catch_Trials  = find (catchtr==1);

    acc           = M1_BehavioralCodes{s}.acc;
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


    TargetSide    = M1_BehavioralCodes{s}.tside;
    TargetSide(TargetSide==2) = -1;
    TargetSide_Used_trials         = TargetSide(Used_Trials);
    TargetSide_Correct_trials      = TargetSide(acc_1);


    % Match_Trials  = zeros(size(TargetSide_Used_trials));
    % Match         = CueSide_Used_trials==TargetSide_Used_trials;
    % idx_match     = find(Match==1);
    % idx_nonmatch  = find(Match==0);
    % Match_Trials    (idx_match)= 1;
    % Match_Trials    (idx_nonmatch)= -1;


    Match_Trials  = zeros(size(An_AC_Sess,1),1);
    Match         = CueLocation==TargetSide;
    idx_match     = find(Match==1);
    idx_nonmatch  = find(Match==0);
    Match_Trials    (idx_match)= 1;
    Match_Trials    (idx_nonmatch)= -1;


    Match_and_Acc=find(sum([Match_Trials acc],2)==2);
    ResponseLabels=zeros(size(Match_Trials))-1 ;
    ResponseLabels(Match_and_Acc)=1;

    % 
    % Match_Corr_Trials  = zeros(size(TargetSide_Correct_trials));
    % Match_Corr         = CueSide_Correct_trials==TargetSide_Correct_trials;
    % idx_Match_Corr     = find(Match_Corr==1);
    % idx_nonMatch_Corr  = find(Match_Corr==0);
    % Match_Corr_Trials    (idx_Match_Corr)= 1;
    % Match_Corr_Trials    (idx_nonMatch_Corr)= -1;


  
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

    PermData_AC=permute(An_AC_Sess,[3 1 2]);
    PermData_PFC=permute(An_PFC_Sess,[3 1 2]);
    % for Del =1:SelectedDelays

    for Del = 4
        if SelectedDelays==1
          Delay_s = 0 ;  
        elseif SelectedDelays==4
        Delay_s = Del ; %%%%  Del
        end
        % Delay = Del ;

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

            Accuracy_AC=[];
            Accuracy_AC=PerformDecodingTrials_AuditoryTask(PermData_AC,Delay_1,Match_Del_1);
            Acc_AC(s,:,Del)=Accuracy_AC;

            Accuracy_PFC=[];
            Accuracy_PFC=PerformDecodingTrials_AuditoryTask(PermData_PFC,Delay_1,Match_Del_1);
            Acc_PFC(s,:,Del)=Accuracy_PFC;

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

            Accuracy_AC=[];
            Accuracy_AC=PerformDecodingTrials_AuditoryTask(PermData_AC,Delay_2,Match_Del_2);
            Acc_AC(s,:,Del)=Accuracy_AC;

            Accuracy_PFC=[];
            Accuracy_PFC=PerformDecodingTrials_AuditoryTask(PermData_PFC,Delay_2,Match_Del_2);
            Acc_PFC(s,:,Del)=Accuracy_PFC;

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

            Accuracy_AC=[];
            Accuracy_AC=PerformDecodingTrials_AuditoryTask(PermData_AC,Delay_3,Match_Del_3);
            Acc_AC(s,:,Del)=Accuracy_AC;

            Accuracy_PFC=[];
            Accuracy_PFC=PerformDecodingTrials_AuditoryTask(PermData_PFC,Delay_3,Match_Del_3);
            Acc_PFC(s,:,Del)=Accuracy_PFC;

        elseif Delay_s  ==4

            PermData_AC=permute(An_AC_Sess,[3 1 2]);
            Accuracy_AC=[];
            Accuracy_AC=PerformDecodingTrials_AuditoryTask(PermData_AC,Delay_All,Match_Trials(Delay_All));
            Acc_AC(s,:,Del)=Accuracy_AC;

            PermData_PFC=permute(An_PFC_Sess,[3 1 2]);
            Accuracy_PFC=[];
            Accuracy_PFC=PerformDecodingTrials_AuditoryTask(PermData_PFC,Delay_All,Match_Trials(Delay_All));
            Acc_PFC(s,:,Del)=Accuracy_PFC;

        end

    end
end


%%


% Times_fig           = -950:25:3550;% for cue alignment
Times_fig           =  -2100:25:1700;% for Tg alignment

options.x_axis      =  Times_fig ;
options.smooth      =  5;
Colors              =  [ [250 160 10]/255; [250 70 10]/255; [4 170 230]/255; ...
    [ 10 90 170]/255;   [60 60 60]/255;  [20 20 20]/255 ];
options.color_PFC   =  Colors(1,:);
options.color_AC    =  Colors(3,:);




d=4
figure('Units', 'pixels', 'Position', [500 500 250 200 ]);

plot_areaerrorbarPFC(Acc_PFC(:,:,d),options); hold on
plot_areaerrorbarAC( Acc_AC(:,:,d),options); hold on
xlim([-500 1000])
ylim([.4 .8])

line([ 0 0 ], [ylim ],'Linestyle','--','Color',[0 0 0 ]/255); hold on









