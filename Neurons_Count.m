clearvars
cd('C:\Users\giarroccof2\OneDrive - National Institutes of Health\Franco\Corrie Data\Task')

load n_neurons_M2
load n_neurons_M1



M=[ n_neurons_M1 n_neurons_M2];


figure('Units', 'pixels', 'Position', [500 500 800 250]);
b = bar(M', 'stacked');
b(1).FaceColor = [4 170 230]/255;
b(2).FaceColor = [250 160 10]/255;
legend({'A1', 'PFC'}, 'Location', 'best','Box','off');
xlabel('Session');
ylabel('Number of Neurons');