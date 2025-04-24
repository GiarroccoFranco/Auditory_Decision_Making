
function Accuracy=PerformDecoding_AuditoryTask(DecData,Trials,Label)

for  bin= 1:size(DecData,3)

    Dat=DecData(:,Trials,bin);
    if ~isempty(Dat)

        DataC=Dat';
        SVMModel = fitcsvm(DataC,Label,"Prior","uniform",'KernelFunction','linear',KFold=10);

        Accuracy(bin) = 1-kfoldLoss(SVMModel);

    elseif isempty(Dat)
        Accuracy(bin) =NaN;
    end

end

end