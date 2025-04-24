
function plot_areaerrorbarPFC(data, options)


    options.alpha      = 0.1;
    options.line_width = 1;
    options.error      = 'sem';
if(isfield(options,'x_axis')==0), options.x_axis = 1:size(data,2); end
options.x_axis = options.x_axis(:);

sm=zeros(size(data));

   sm=smoothdata(data,2,'gaussian',options.s );

data=sm;

data_mean = nanmean(data,1);
for l=1:size(data,2)
    bin=[];
    bin=data(:,l);
    idx = isnan(bin);
    bin = bin(~idx);
    data_std(l)=std(bin);
end


error = (data_std./sqrt(size(data,1)));

x_vector = [options.x_axis', fliplr(options.x_axis')];
patch = fill(x_vector, [data_mean+error,fliplr(data_mean-error)], options.color_PFC);
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', options.alpha);
hold on;
plot(options.x_axis, data_mean, 'color', options.color_PFC, ...
    'LineWidth', options.line_width);

end