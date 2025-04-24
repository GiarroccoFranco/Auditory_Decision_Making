%%% M1


load  Micro_Diff_1lag_Peak_3
load  Micro_MnM_PeakTime_1lag
load  Micro_Latency_MnM_1lag

load  Micro_Diff_2lag_Peak_3
load  Micro_MnM_PeakTime_2lag
load  Micro_Latency_MnM_2lag

load  Micro_Diff_3lag_Peak_3
load  Micro_MnM_PeakTime_3lag
load  Micro_Latency_MnM_3lag

load  Micro_Diff_4lag_Peak_3
load  Micro_MnM_PeakTime_4lag
load  Micro_Latency_MnM_4lag

%%% M2

load  Mega_Diff_1lag_Peak_3
load  Mega_MnM_PeakTime_1lag
load  Mega_Latency_MnM_1lag

load  Mega_Diff_2lag_Peak_3
load  Mega_MnM_PeakTime_2lag
load  Mega_Latency_MnM_2lag

load  Mega_Diff_3lag_Peak_3
load  Mega_MnM_PeakTime_3lag
load  Mega_Latency_MnM_3lag

load  Mega_Diff_4lag_Peak_3
load  Mega_MnM_PeakTime_4lag
load  Mega_Latency_MnM_4lag

AcPeaks=[   vertcat(Micro_Diff_1lag_Peak_3(:,1), Mega_Diff_1lag_Peak_3(:,1)) , ...
            vertcat(Micro_Diff_2lag_Peak_3(:,1),      Mega_Diff_2lag_Peak_3(:,1))      , ...
            vertcat(Micro_Diff_3lag_Peak_3(:,1), Mega_Diff_3lag_Peak_3(:,1)) , ...
            vertcat(Micro_Diff_4lag_Peak_3(:,1), Mega_Diff_4lag_Peak_3(:,1))];


PFCPeaks=[   vertcat(Micro_Diff_1lag_Peak_3(:,2), Mega_Diff_1lag_Peak_3(:,2)) , ...
            vertcat(Micro_Diff_2lag_Peak_3(:,2),      Mega_Diff_2lag_Peak_3(:,2))      , ...
            vertcat(Micro_Diff_3lag_Peak_3(:,2), Mega_Diff_3lag_Peak_3(:,2)) , ...
            vertcat(Micro_Diff_4lag_Peak_3(:,2), Mega_Diff_4lag_Peak_3(:,2))];


AcPeaks = abs(AcPeaks);
PFCPeaks = abs(PFCPeaks);


%%

%%
% % % % % For Peak

Mega_Peaks= { [Mega_MnM_PeakTime_1lag{1} Mega_MnM_PeakTime_1lag{2}]  ,...
              [Mega_MnM_PeakTime_2lag{1}      Mega_MnM_PeakTime_2lag{2}]  ,...
              [Mega_MnM_PeakTime_3lag{1} Mega_MnM_PeakTime_3lag{2}]  ,...
              [Mega_MnM_PeakTime_4lag{1} Mega_MnM_PeakTime_4lag{2}]} ;

Micro_Peaks= { [Micro_MnM_PeakTime_1lag{1} Micro_MnM_PeakTime_1lag{2}]  ,...
              [Micro_MnM_PeakTime_2lag{1}      Micro_MnM_PeakTime_2lag{2}]  ,...
              [Micro_MnM_PeakTime_3lag{1} Micro_MnM_PeakTime_3lag{2}]  ,...
              [Micro_MnM_PeakTime_4lag{1} Micro_MnM_PeakTime_4lag{2}]} ;

IT_Peaks = {Micro_Peaks, Mega_Peaks};

for MK=1:2
Latenc=IT_Peaks{MK};

for lags= 1: size (Latenc,2)

    lag_latenc=[]; lag_latenc=Latenc{lags};

    [~,p_lag_Peak(lags,MK)] = ttest(lag_latenc(:,1),lag_latenc(:,2));

    Mean_Peak(lags,:,MK)=mean(lag_latenc);
    SEM_Peak(lags,:,MK)=std(lag_latenc)/sqrt(size(lag_latenc,1));

end


PFC_Latency = [Latenc{1}(:,1) Latenc{2}(:,1) Latenc{3}(:,1) Latenc{4}(:,1)];
A1_Latency = [Latenc{1}(:,2) Latenc{2}(:,2) Latenc{3}(:,2) Latenc{4}(:,2)];


Average_peak_and_SEM_PFC(1,MK)=mean(PFC_Latency,'all');
Average_peak_and_SEM_PFC(2,MK)=std(mean(PFC_Latency,2))/sqrt(size(PFC_Latency,1));


Average_peak_and_SEM_A1(1,MK)=mean(A1_Latency,'all');
Average_peak_and_SEM_A1(2,MK)=std(mean(A1_Latency,2))/sqrt(size(A1_Latency,1));


Y_latency_4_Anova=cat(2,PFC_Latency,A1_Latency) ;
Group_Area1=zeros(size(PFC_Latency));
Group_Area2=ones(size(A1_Latency)); 
Areas_Group=cat(2,Group_Area1,Group_Area2) ;

Prediction_Bins=ones(size(Areas_Group));
Prediction_Bins(:,[2 6])=2;
Prediction_Bins(:,[3 7])=3;
Prediction_Bins(:,[4 8])=4;
anovan(Y_latency_4_Anova(:),[Areas_Group(:) Prediction_Bins(:)],'model','full');


Mean_Peak_RT_Diff_PFC(MK)=mean(RTs_UsedTrials{MK}*1000-Average_peak_and_SEM_PFC(1,MK))
Sem_Peak_RT_Diff_PFC(MK)=std(RTs_UsedTrials{MK}*1000-Average_peak_and_SEM_PFC(1,MK))/numel(RTs_UsedTrials{MK});

Mean_Peak_RT_Diff_A1(MK)=mean(RTs_UsedTrials{MK}*1000-Average_peak_and_SEM_A1(1,MK))
Sem_Peak_RT_Diff_A1(MK)=std(RTs_UsedTrials{MK}*1000-Average_peak_and_SEM_A1(1,MK))/numel(RTs_UsedTrials{MK});

end

figure('Units', 'pixels', 'Position', [500 500 700 300]);
for MK=1:2

subplot(1,2,MK)
errorbar(1:4,Mean_Peak(:,1,MK),SEM_Peak(:,1,MK),'Color',[250 160 10]/255 ,"LineStyle","-"); hold on
errorbar(1:4,Mean_Peak(:,2,MK),SEM_Peak(:,2,MK),'Color',[4 170 230]/255 ,"LineStyle","-"); hold on



ylim([400 650])
xlim([0 5])
end




%%
%%%%% For latency

Mega_latencies= { Mega_Latency_MnM_1lag, Mega_Latency_MnM_2lag , Mega_Latency_MnM_3lag, Mega_Latency_MnM_4lag} ;

Micro_latencies= { Micro_Latency_MnM_1lag, Micro_Latency_MnM_2lag , Micro_Latency_MnM_3lag, Micro_Latency_MnM_4lag} ;

IT_Latencies = { Micro_latencies, Mega_latencies };


for MK=1:2

Latenc=IT_Latencies{MK};
AllLags_Mean=[];
for lags= 1: size (Latenc,2)

    lag_latenc=[]; lag_latenc=Latenc{lags};

    [~,p_lag_latency(lags,MK)] = ttest2(lag_latenc(:,1),lag_latenc(:,2));
    for ar=1:2
    AllLags_Mean  (:,lags,ar)=lag_latenc(:,ar);


    end
    Mean_Latency(lags,:,MK)=mean(lag_latenc);
    SEM_Latency(lags,:,MK)=std(lag_latenc)/sqrt(size(lag_latenc,1));

end


mean_Lat=mean(AllLags_Mean,2);
mean_Latency=mean(mean_Lat);  SEM_Lat=std(mean_Lat)/sqrt(size(mean_Lat,1))
[l1,l2,l3,l4]=ttest(mean_Lat(:,:,1),mean_Lat(:,:,2))
PFC_Latency = [Latenc{1}(:,1) Latenc{2}(:,1) Latenc{3}(:,1) Latenc{4}(:,1)];
A1_Latency = [Latenc{1}(:,2) Latenc{2}(:,2) Latenc{3}(:,2) Latenc{4}(:,2)];


Mean_Latency_RT_Diff(:,:,MK)=mean(RTs_UsedTrials{MK}*1000-mean_Lat)
Sem_Latency_RT_Diff(:,:,MK)=std(RTs_UsedTrials{MK}*1000-mean_Lat)/numel(RTs_UsedTrials{MK});


Average_peak_and_SEM_PFC(1,MK)=mean(PFC_Latency,'all');
Average_peak_and_SEM_PFC(2,MK)=std(mean(PFC_Latency,2))/sqrt(size(PFC_Latency,1));

Average_peak_and_SEM_A1(1,MK)=mean(A1_Latency,'all');
Average_peak_and_SEM_A1(2,MK)=std(mean(A1_Latency,2))/sqrt(size(A1_Latency,1));


Y_latency_4_Anova=cat(2,PFC_Latency,A1_Latency) ;
Group_Area1=zeros(size(PFC_Latency));
Group_Area2=ones(size(A1_Latency)); 
Areas_Group=cat(2,Group_Area1,Group_Area2) ;

Prediction_Bins=ones(size(Areas_Group));
Prediction_Bins(:,[2 6])=2;
Prediction_Bins(:,[3 7])=3;
Prediction_Bins(:,[4 8])=4;


Sessions=1:size(Prediction_Bins,1);
groupSession=[];
for b=1:size(Prediction_Bins,2)
groupSession(:,b)=Sessions;
end


[p1, tbl1, stats1] = anovan(Y_latency_4_Anova(:),[Areas_Group(:) Prediction_Bins(:) ]);

end

figure('Units', 'pixels', 'Position', [500 500 700 300]);

for MK=1:2

subplot(1,2,MK)
errorbar(1:4,Mean_Latency(:,1,MK),SEM_Latency(:,1,MK),'Color',[250 160 10]/255 ,"LineStyle","-"); hold on
errorbar(1:4,Mean_Latency(:,2,MK),SEM_Latency(:,2,MK),'Color',[4 170 230]/255 ,"LineStyle","-"); hold on

xlim([0 5])
if MK==1
    ylim([250 600])

elseif MK==2
ylim([250 650])
end

end

