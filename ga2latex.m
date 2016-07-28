function [] = ga2latex(SubStruct,dir_ref,name)

cleanName = strrep(name, '.', '');
cleanName = strrep(cleanName, ' ', '_');

filename = ['LaTexTables_Global_',cleanName,'.txt'];
filepath = [dir_ref,'/',filename];

lines = 30;
tobuild = cell(lines,1);

structureNames = {'KolD_Pri';'CvM_Pri';'Kuiper_Pri';'Watson_Pri';'AD_Pri'};
labels = {'Priority: Kolmogorov-D','Priority: Cramer-von-Mises','Priority: Kuiper','Priority: Watson','Priority: Anderson-Darling'};

paraCount = zeros(1,length(structureNames));

distNames = cell(1,length(structureNames));

parfor i=1:length(structureNames)
    thisName = structureNames{i};
    thisDist = SubStruct.(thisName).Distribution;
    thisPara = length(thisDist)-1;
    paraCount(i) = thisPara;
    distNames{i} = thisDist.Type;
end

upperParaMax = max(paraCount(1:3));
lowerParaMax = max(paraCount(4:end));

for i=1:length(structureNames)
    thisName = structureNames{i};
    thisDist = SubStruct.(thisName).Distribution;
    if i>3
        height = lowerParaMax;
    else
        height = upperParaMax;
    end
    para.(thisName) = parameterLaTeX(thisDist,height);
end

if mod(length(structureNames),2)
    blankstruc = struct('Type',[]);
    para.blank = parameterLaTeX(blankstruc,lowerParaMax);
    paraCount = [paraCount 1];
end

clineUpper = clineLaTex(paraCount(1:3),upperParaMax);
clineLower = clineLaTex(paraCount(4:end),lowerParaMax);

%Construct LaTeX
tobuild{001} = '\begin{tabular}{|c|c||c|c||c|c||c|c|} \hline';
tobuild{002} = ['\multicolumn{2}{|c||}{',name,'} & \multicolumn{2}{c||}{\textbf{',labels{1},'}} & \multicolumn{2}{c||}{\textbf{',labels{2},'}} & \multicolumn{2}{c||}{\textbf{',labels{3},'}} \\ \hline'];
tobuild{003} = ['\multicolumn{2}{|c||}{\textbf{Best Distribution}} & \multicolumn{2}{c||}{\textbf{',distNames{1},'}} & \multicolumn{2}{c||}{\textbf{',distNames{2},'}} & \multicolumn{2}{c||}{\textbf{',distNames{3},'}} \\ \hline'];
if upperParaMax == 3
    tobuild{004} = ['\multicolumn{2}{|c||}{\multirow{3}{*}{\textbf{Parameters}}} & ',para.KolD_Pri{1},' & ',para.CvM_Pri{1},' & ',para.Kuiper_Pri{1},'\\ ',clineUpper{1}];
    tobuild{005} = ['\multicolumn{2}{|c||}{} & ',para.KolD_Pri{2},' & ',para.CvM_Pri{2},' & ',para.Kuiper_Pri{2},'\\ ',clineUpper{2}];
    tobuild{006} = ['\multicolumn{2}{|c||}{} & ',para.KolD_Pri{3},' & ',para.CvM_Pri{3},' & ',para.Kuiper_Pri{3},'\\ ',clineUpper{3}];
elseif upperParaMax == 2
    tobuild{004} = ['\multicolumn{2}{|c||}{\multirow{2}{*}{\textbf{Parameters}}} & ',para.KolD_Pri{1},' & ',para.CvM_Pri{1},' & ',para.Kuiper_Pri{1},'\\ ',clineUpper{1}];
    tobuild{005} = ['\multicolumn{2}{|c||}{} & ',para.KolD_Pri{2},' & ',para.CvM_Pri{2},' & ',para.Kuiper_Pri{2},'\\ ',clineUpper{2}];
    tobuild{006} = '%This line has been intentionally left blank';    
elseif upperParamax == 1
    tobuild{004} = ['\multicolumn{2}{|c||}{\textbf{Parameters}} & ',para.KolD_Pri(1),' & ',para.CvM_Pri(1),' & ',para.Kuiper_Pri(1),'\\ ',clineUpper{1}];
    tobuild{005} = '%This line has been intentionally left blank';
    tobuild{006} = '%This line has been intentionally left blank';
end
tobuild{007} = ['\parbox[t]{2mm}{\multirow{5}{*}{\rotatebox[origin=c]{90}{Statistics}}} & Kol D & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.Statistics.Kolmogorov_D),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.Statistics.Kolmogorov_D),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.Statistics.Kolmogorov_D),'$} \\ \cline{2-8}'];
tobuild{008} = ['& CvM & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.Statistics.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.Statistics.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.Statistics.Cramer_von_Mises),'$} \\ \cline{2-8}'];
tobuild{009} = ['& Kuiper & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.Statistics.Kuiper),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.Statistics.Kuiper),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.Statistics.Kuiper),'$} \\ \cline{2-8}'];
tobuild{010} = ['& Watson & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.Statistics.Watson),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.Statistics.Watson),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.Statistics.Watson),'$} \\ \cline{2-8}'];
tobuild{011} = ['& And-Dar & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.Statistics.Anderson_Darling),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.Statistics.Anderson_Darling),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.Statistics.Anderson_Darling),'$} \\ \hline'];
tobuild{012} = ['\parbox[t]{2mm}{\multirow{5}{*}{\rotatebox[origin=c]{90}{p-Values}}} & Kol D & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.pValues.Kolmogorov_D),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.pValues.Kolmogorov_D),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.pValues.Kolmogorov_D),'$} \\ \cline{2-8}'];
tobuild{013} = ['& CvM & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.pValues.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.pValues.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.pValues.Cramer_von_Mises),'$} \\ \cline{2-8}'];
tobuild{014} = ['& Kuiper & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.pValues.Kuiper),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.pValues.Kuiper),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.pValues.Kuiper),'$} \\ \cline{2-8}'];
tobuild{015} = ['& Watson & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.pValues.Watson),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.pValues.Watson),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.pValues.Watson),'$} \\ \cline{2-8}'];
tobuild{016} = ['& And-Dar & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.pValues.Anderson_Darling),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.pValues.Anderson_Darling),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.pValues.Anderson_Darling),'$} \\ \hline\hline'];

tobuild{017} = ['\multicolumn{2}{|c||}{\textbf{Best Distribution}} & \multicolumn{2}{c||}{\textbf{',distNames{4},'}} & \multicolumn{2}{c||}{\textbf{',distNames{5},'}} & \multicolumn{2}{c||}{\textbf{}} \\ \hline'];
if lowerParaMax == 3
    tobuild{018} = ['\multicolumn{2}{|c||}{\multirow{3}{*}{\textbf{Parameters}}} & ',para.Watson_Pri{1},' & ',para.AD_Pri{1},' & ',para.Blank{1},'\\ ',clineLower{1}];
    tobuild{019} = ['\multicolumn{2}{|c||}{} & ',para.Watson_Pri{2},' & ',para.AD_Pri{2},' & ',para.Blank{2},'\\ ',clineLower{2}];
    tobuild{020} = ['\multicolumn{2}{|c||}{} & ',para.Watson_Pri{3},' & ',para.AD_Pri{3},' & ',para.Blank{3},'\\ ',clineLower{3}];
elseif lowerParaMax == 2
    tobuild{018} = ['\multicolumn{2}{|c||}{\multirow{2}{*}{\textbf{Parameters}}} & ',para.Watson_Pri{1},' & ',para.AD_Pri{1},' & ',para.Blank{1},'\\ ',clineLower{1}];
    tobuild{019} = ['\multicolumn{2}{|c||}{} & ',para.Watson_Pri{2},' & ',para.AD_Pri{2},' & ',para.Blank{2},'\\ ',clineLower{2}];
    tobuild{020} = '%This line has been intentionally left blank';    
elseif lowerParaMax == 1
    tobuild{018} = ['\multicolumn{2}{|c||}{\textbf{Parameters}} & ',para.Watson_Pri(1),' & ',para.AD_Pri(1),' & ',para.Blank(1),'\\ ',clineLower{1}];
    tobuild{019} = '%This line has been intentionally left blank';
    tobuild{020} = '%This line has been intentionally left blank';
end
tobuild{021} = ['\parbox[t]{2mm}{\multirow{5}{*}{\rotatebox[origin=c]{90}{Statistics}}} & Kol D & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.Statistics.Kolmogorov_D),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.Statistics.Kolmogorov_D),'$} & \multicolumn{2}{c||}{} \\ \cline{2-8}'];
tobuild{022} = ['& CvM & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.Statistics.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.Statistics.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{} \\ \cline{2-8}'];
tobuild{023} = ['& Kuiper & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.Statistics.Kuiper),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.Statistics.Kuiper),'$} & \multicolumn{2}{c||}{} \\ \cline{2-8}'];
tobuild{024} = ['& Watson & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.Statistics.Watson),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.Statistics.Watson),'$} & \multicolumn{2}{c||}{} \\ \cline{2-8}'];
tobuild{025} = ['& And-Dar & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.Statistics.Anderson_Darling),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.Statistics.Anderson_Darling),'$} & \multicolumn{2}{c||}{} \\ \hline'];
tobuild{026} = ['\parbox[t]{2mm}{\multirow{5}{*}{\rotatebox[origin=c]{90}{p-Values}}} & Kol D & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.pValues.Kolmogorov_D),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.pValues.Kolmogorov_D),'$} & \multicolumn{2}{c||}{} \\ \cline{2-8}'];
tobuild{027} = ['& CvM & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.pValues.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.pValues.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{} \\ \cline{2-8}'];
tobuild{028} = ['& Kuiper & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.pValues.Kuiper),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.pValues.Kuiper),'$} & \multicolumn{2}{c||}{} \\ \cline{2-8}'];
tobuild{029} = ['& Watson & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.pValues.Watson),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.pValues.Watson),'$} & \multicolumn{2}{c||}{} \\ \cline{2-8}'];
tobuild{030} = ['& And-Dar & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.pValues.Anderson_Darling),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.pValues.Anderson_Darling),'$} & \multicolumn{2}{c||}{} \\ \hline'];

fileID = fopen(filepath,'w');
fprintf(fileID,'%s\r\n',tobuild{:});
fclose(fileID);

end