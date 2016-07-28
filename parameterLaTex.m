function [insert] = parameterLaTex(distStruct,height)

insert = cell(height,1);

if height == 1
    if strcmp(distStruct.Type,'Exponential')
        insert{1} = ['Scale & $',num2matlabstr(distStruct.Scale),'$'];
    elseif strcmp(distStruct.Type,'Rayleigh')
        insert{1} = ['Shape & $',num2matlabstr(distStruct.Shape),'$'];
    else
        insert{1} = ' & ';
    end
elseif height == 2
    if strcmp(distStruct.Type,'Exponential')
        insert{1} = ['\multirow{2}{*}{Scale} & \multirow{2}{*}{$',num2matlabstr(distStruct.Scale),'$}'];
        insert{2} = ' & ';
    elseif strcmp(distStruct.Type,'Gamma')
        insert{1} = ['Shape & $',num2matlabstr(distStruct.Shape),'$'];
        insert{2} = ['Scale & $',num2matlabstr(distStruct.Scale),'$'];
    elseif strcmp(distStruct.Type,'Rayleigh')
        insert{1} = ['\multirow{2}{*}{Scale} & \multirow{2}{*}{$',num2matlabstr(distStruct.Scale),'$}'];
        insert{2} = ' & ';
    elseif strcmp(distStruct.Type,'Log-Normal')
        insert{1} = ['Shape & $',num2matlabstr(distStruct.Shape),'$'];
        insert{2} = ['Scale & $',num2matlabstr(distStruct.Scale),'$'];
    elseif strcmp(distStruct.Type,'Mittag-Leffler')
        insert{1} = ['Stability & $',num2matlabstr(distStruct.Stability),'$'];
        insert{2} = ['Scale & $',num2matlabstr(distStruct.Scale),'$'];
    elseif strcmp(distStruct.Type,'Weibull')
        insert{1} = ['Scale & $',num2matlabstr(distStruct.Scale),'$'];
        insert{2} = ['Shape & $',num2matlabstr(distStruct.Shape),'$'];
    else
        insert{1} = ' & ';
        insert{2} = ' & ';
    end
elseif height == 3
    if strcmp(distStruct.Type,'Exponential')
        insert{1} = ['\multirow{3}{*}{Scale} & \multirow{3}{*}{$',num2matlabstr(distStruct.Scale),'$}'];
        insert{2} = ' & ';
        insert{3} = ' & ';
    elseif strcmp(distStruct.Type,'Gamma')
        insert{1} = ['Shape & $',num2matlabstr(distStruct.Shape),'$'];
        insert{2} = ['\multirow{2}{*}{Scale} & \multirow{2}{*}{$',num2matlabstr(distStruct.Scale),'$}'];
        insert{3} = ' & ';
    elseif strcmp(distStruct.Type,'Rayleigh')
        insert{1} = ['\multirow{3}{*}{Scale} & \multirow{3}{*}{$',num2matlabstr(distStruct.Scale),'$}'];
        insert{2} = ' & ';
        insert{3} = ' & ';
    elseif strcmp(distStruct.Type,'Log-Normal')
        insert{1} = ['Shape & $',num2matlabstr(distStruct.Shape),'$'];
        insert{2} = ['\multirow{2}{*}{Scale} & \multirow{2}{*}{$',num2matlabstr(distStruct.Scale),'$}'];
        insert{3} = ' & ';
    elseif strcmp(distStruct.Type,'Mittag-Leffler')
        insert{1} = ['Stability & $',num2matlabstr(distStruct.Stability),'$'];
        insert{2} = ['\multirow{2}{*}{Scale} & \multirow{2}{*}{$',num2matlabstr(distStruct.Scale),'$}'];
        insert{3} = ' & ';
    elseif strcmp(distStruct.Type,'Generalised Pareto')
        insert{1} = ['Shape & $',num2matlabstr(distStruct.Shape),'$'];
        insert{2} = ['Scale & $',num2matlabstr(distStruct.Scale),'$'];
        insert{3} = ['Location & $',num2matlabstr(distStruct.Location),'$'];
    elseif strcmp(distStruct.Type,'Weibull')
        insert{1} = ['Scale & $',num2matlabstr(distStruct.Scale),'$'];
        insert{2} = ['\multirow{2}{*}{Shape} & \multirow{2}{*}{$',num2matlabstr(distStruct.Shape),'$}'];
        insert{3} = ' & ';
    else
        insert{1} = ' & ';
        insert{2} = ' & ';
        insert{3} = ' & ';
    end
end