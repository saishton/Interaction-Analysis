function [] = global_plot(data,bestFits,dir_ref,property_title,graph_title,cutExtreme,Ymin)

fieldNames = fieldnames(data);
numData = length(fieldNames);
structureNames = {'KolD_Pri';'CvM_Pri';'Kuiper_Pri';'Watson_Pri';'AD_Pri';'KL_Pri';'JS_Pri';};
labels = {'Data','Priority: Kolmogorov-D','Priority: Cramer-von-Mises','Priority: Kuiper','Priority: Watson','Priority: Anderson-Darling','Priority: Kullback-Leibler','Priority: Jensen-Shannon'};

splices = 1000;
Xmax = 0;
h = zeros(1,numData+length(structureNames));

fig = figure();
hold on

for p=1:numData
    thisname = fieldNames{p};
    thisdata = data.(thisname);
    [F,X] = ecdf(thisdata);
    ccdf = 1-F;
    X = X(2:end-cutExtreme);
    ccdf = ccdf(2:end-cutExtreme);
    Xmax = max(Xmax,X(end));
    h(p) = plot(X,ccdf,'o');
end

grid = linspace(0,Xmax,splices);
for q=1:length(structureNames)
    thisName = structureNames{q};
    thisDistribution = bestFits.(thisName).Distribution;
    if strcmp(thisDistribution.Type,'Exponential')
        ccdf = expcdf(grid,thisDistribution.Scale,'upper');
    elseif strcmp(thisDistribution.Type,'Gamma')
        ccdf = gamcdf(grid,thisDistribution.Shape,thisDistribution.Scale,'upper');
    elseif strcmp(thisDistribution.Type,'Rayleigh')
        ccdf = raylcdf(grid,thisDistribution.Shape,'upper');
    elseif strcmp(thisDistribution.Type,'Log-Normal')
        ccdf = logncdf(grid,thisDistribution.Shape,thisDistribution.Scale,'upper');
    elseif strcmp(thisDistribution.Type,'Mittag-Leffler')
        ccdf = mlf(thisDistribution.Stability,1,-thisDistribution.Scale*grid.^thisDistribution.Stability,6);
    elseif strcmp(thisDistribution.Type,'Generalised Pareto')
        ccdf = gpcdf(grid,thisDistribution.Shape,thisDistribution.Scale,thisDistribution.Location,'upper');
    elseif strcmp(thisDistribution.Type,'Weibull')
        ccdf = wblcdf(grid,thisDistribution.Scale,thisDistribution.Shape,'upper');
    end
    h(numData+q) = plot(grid,ccdf);
end

set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel(graph_title);
ylabel('CCDF');
axis([-inf,inf,Ymin,1E0]);
h = [h(1) h(numData+1:end)];
legend(h,labels,'Location','southwest');
imagefilename = [dir_ref,'/',property_title,'.png'];
print(imagefilename,'-dpng')
close(fig);