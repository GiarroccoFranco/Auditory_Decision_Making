load Mega_Cue_PeakValue_1lag
load Mega_Latency_Cue_1lag
load Mega_Cue_PeakTime_1lag

load Mega_Cue_PeakValue_2lag
load Mega_Latency_Cue_2lag
load Mega_Cue_PeakTime_2lag

load Mega_Cue_PeakValue_3lag
load Mega_Latency_Cue_3lag
load Mega_Cue_PeakTime_3lag

load Mega_Cue_PeakValue_4lag
load Mega_Latency_Cue_4lag
load Mega_Cue_PeakTime_4lag

%-----------------------------

load Micro_Cue_PeakValue_1lag
load Micro_Latency_Cue_1lag
load Micro_Cue_PeakTime_1lag

load Micro_Cue_PeakValue_2lag
load Micro_Latency_Cue_2lag
load Micro_Cue_PeakTime_2lag

load Micro_Cue_PeakValue_3lag
load Micro_Latency_Cue_3lag
load Micro_Cue_PeakTime_3lag

load Micro_Cue_PeakValue_4lag
load Micro_Latency_Cue_4lag
load Micro_Cue_PeakTime_4lag




%% for Peaks

Micro_CuePeak =  [ Micro_Cue_PeakTime_1lag{2}  Micro_Cue_PeakTime_2lag{2} Micro_Cue_PeakTime_3lag{2} Micro_Cue_PeakTime_4lag{2}];
Mega_CuePeak =  [ Mega_Cue_PeakTime_1lag{2}  Mega_Cue_PeakTime_2lag{2} Mega_Cue_PeakTime_3lag{2} Mega_Cue_PeakTime_4lag{2}];

Cue_Peak= { Micro_CuePeak, Mega_CuePeak} ;

figure('Units', 'pixels', 'Position', [500 500 600 250]);

for MK=1:2
MK_Peak=Cue_Peak{MK};

Mean_Latency=mean(MK_Peak);
Sem_Latency= std(MK_Peak)/sqrt(size(MK_Peak(:,1),1));
subplot(1,2,MK)
errorbar([1:4], Mean_Latency,Sem_Latency ); hold on
ylim([50 200])
xlim([0 5])

Average_peak_and_SEM(1,MK)=mean(MK_Peak,'all');
Average_peak_and_SEM(2,MK)=std(mean(MK_Peak,2))/sqrt(size(MK_Peak,1));
Anova_Group=zeros(size(MK_Peak)); Anova_Group(:,1)=1;Anova_Group(:,2)=2;Anova_Group(:,3)=3;Anova_Group(:,4)=4;
Sessions=1:size(MK_Peak,1);
groupSession=[];
for b=1:size(MK_Peak,2)
groupSession(:,b)=Sessions;
end

anovan(MK_Peak(:),{Anova_Group(:), groupSession(:) },'random',2);

end

%% for latency

Micro_CueLatency = [ Micro_Latency_Cue_1lag(:,2)  Micro_Latency_Cue_2lag(:,2) Micro_Latency_Cue_3lag(:,2) Micro_Latency_Cue_4lag(:,2)];
Mega_CueLatency = [ Mega_Latency_Cue_1lag(:,2)  Mega_Latency_Cue_2lag(:,2) Mega_Latency_Cue_3lag(:,2) Mega_Latency_Cue_4lag(:,2)];

Cue_Latency= { Micro_CueLatency, Mega_CueLatency} ;
figure('Units', 'pixels', 'Position', [500 500 600 250]);

for MK=1:2
MK_Latency=Cue_Latency{MK};

Mean_Latency=mean(MK_Latency);
Sem_Latency= std(MK_Latency)/sqrt(size(MK_Latency(:,1),1));
subplot(1,2,MK)
errorbar([1:4], Mean_Latency,Sem_Latency ); hold on
ylim([50 200])
xlim([0 5])


Sessions=1:size(MK_Peak,1);
groupSession=[];
for b=1:size(MK_Peak,2)
groupSession(:,b)=Sessions;
end

Average_latency_and_SEM(1,MK)=mean(MK_Latency,'all');
Average_latency_and_SEM(2,MK)=std(mean(MK_Latency,2))/sqrt(size(MK_Latency,1));
Anova_Group=zeros(size(MK_Latency)); Anova_Group(:,1)=1;Anova_Group(:,2)=2;Anova_Group(:,3)=3;Anova_Group(:,4)=4;
Sessions=1:size(MK_Latency,1);
groupSession=[];

for b=1:size(MK_Latency,2)
groupSession(:,b)=Sessions;
end

Groups={ Anova_Group(:) groupSession(:)   };
Varnames={'# of bin','session'  };


anovan(MK_Latency(:),Groups    ,'random',2,'varnames',Varnames, model='interaction');

end


   
