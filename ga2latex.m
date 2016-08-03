function [] = ga2latex(SubStruct,dir_ref,name)

cleanName = strrep(name, '.', '');
cleanName = strrep(cleanName, ' ', '_');

filename = ['LaTexTables_Global_',cleanName,'.txt'];
filepath = [dir_ref,'/',filename];

lines = 39;
tobuild = cell(lines,1);

structureNames = {'KolD_Pri';'CvM_Pri';'Kuiper_Pri';'Watson_Pri';'AD_Pri';'KL_Pri';'JS_Pri';};
labels = {'Priority: Kolmogorov-D','Priority: Cramer-von-Mises','Priority: Kuiper','Priority: Watson','Priority: Anderson-Darling','Priority: Kullback-Leibler','Priority: Jensen-Shannon'};

paraCount = zeros(1,length(structureNames));

distNames = cell(1,length(structureNames));

parfor i=1:length(structureNames)
    thisName = structureNames{i};
    thisDist = SubStruct.(thisName).Distribution;
    thisPara = length(thisDist)-1;
    paraCount(i) = thisPara;
    distNames{i} = thisDist.Type;
end

upperParaMax = max(paraCount(1:4));
lowerParaMax = max(paraCount(5:end));

for i=1:length(structureNames)
    thisName = structureNames{i};
    thisDist = SubStruct.(thisName).Distribution;
    if i>4
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

clineUpper = clineLaTex(paraCount(1:4),upperParaMax);
clineLower = clineLaTex(paraCount(5:end),lowerParaMax);

%Construct LaTeX
tobuild{001} = '\begin{tabular}{|c|c||c|c||c|c||c|c||c|c|} \hline';
tobuild{002} = ['\multicolumn{2}{|c||}{',name,'} & \multicolumn{2}{c||}{\textbf{',labels{1},'}} & \multicolumn{2}{c||}{\textbf{',labels{2},'}} & \multicolumn{2}{c||}{\textbf{',labels{3},'}} & \multicolumn{2}{c||}{\textbf{',labels{4},'}} \\ \hline'];
tobuild{003} = ['\multicolumn{2}{|c||}{\textbf{Best Distribution}} & \multicolumn{2}{c||}{\textbf{',distNames{1},'}} & \multicolumn{2}{c||}{\textbf{',distNames{2},'}} & \multicolumn{2}{c||}{\textbf{',distNames{3},'}} & \multicolumn{2}{c||}{\textbf{',distNames{4},'}} \\ \hline'];
if upperParaMax == 3
    tobuild{004} = ['\multicolumn{2}{|c||}{\multirow{3}{*}{\textbf{Parameters}}} & ',para.A_Pri{1},' & ',para.CvM_Pri{1},' & ',para.Kuiper_Pri{1},' & ',para.Watson_Pri{1},'\\ ',clineUpper{1}];
    tobuild{005} = ['\multicolumn{2}{|c||}{} & ',para.KolD_Pri{2},' & ',para.CvM_Pri{2},' & ',para.Kuiper_Pri{2},' & ',para.Watson_Pri{2},'\\ ',clineUpper{2}];
    tobuild{006} = ['\multicolumn{2}{|c||}{} & ',para.KolD_Pri{3},' & ',para.CvM_Pri{3},' & ',para.Kuiper_Pri{3},' & ',para.Watson_Pri{3},'\\ ',clineUpper{3}];
elseif upperParaMax == 2
    tobuild{004} = ['\multicolumn{2}{|c||}{\multirow{2}{*}{\textbf{Parameters}}} & ',para.KolD_Pri{1},' & ',para.CvM_Pri{1},' & ',para.Kuiper_Pri{1},' & ',para.Watson_Pri{1},'\\ ',clineUpper{1}];
    tobuild{005} = ['\multicolumn{2}{|c||}{} & ',para.KolD_Pri{2},' & ',para.CvM_Pri{2},' & ',para.Kuiper_Pri{2},' & ',para.Watson_Pri{2},'\\ ',clineUpper{2}];
    tobuild{006} = '%This line has been intentionally left blank';    
elseif upperParamax == 1
    tobuild{004} = ['\multicolumn{2}{|c||}{\textbf{Parameters}} & ',para.KolD_Pri(1),' & ',para.CvM_Pri(1),' & ',para.Kuiper_Pri{1},' & ',para.Watson_Pri{1},'\\ ',clineUpper{1}];
    tobuild{005} = '%This line has been intentionally left blank';
    tobuild{006} = '%This line has been intentionally left blank';
end
tobuild{007} = ['\parbox[t]{2mm}{\multirow{7}{*}{\rotatebox[origin=c]{90}{Statistics}}} & Kol D & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.Statistics.Kolmogorov_D),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.Statistics.Kolmogorov_D),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.Statistics.Kolmogorov_D),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.Statistics.Kolmogorov_D),'$} \\ \cline{2-10}'];
tobuild{008} = ['& CvM & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.Statistics.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.Statistics.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.Statistics.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.Statistics.Cramer_von_Mises),'$} \\ \cline{2-10}'];
tobuild{009} = ['& Kuiper & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.Statistics.Kuiper),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.Statistics.Kuiper),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.Statistics.Kuiper),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.Statistics.Kuiper),'$} \\ \cline{2-10}'];
tobuild{010} = ['& Watson & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.Statistics.Watson),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.Statistics.Watson),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.Statistics.Watson),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.Statistics.Watson),'$} \\ \cline{2-10}'];
tobuild{011} = ['& And-Dar & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.Statistics.Anderson_Darling),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.Statistics.Anderson_Darling),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.Statistics.Anderson_Darling),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.Statistics.Anderson_Darling),'$} \\ \cline{2-10}'];
tobuild{012} = ['& Kul-Lei & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.Statistics.Kullback_Leibler),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.Statistics.Kullback_Leibler),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.Statistics.Kullback_Leibler),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.Statistics.Kullback_Leibler),'$} \\ \cline{2-10}'];
tobuild{013} = ['& Jen-Sha & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.Statistics.Jensen_Shannon),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.Statistics.Jensen_Shannon),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.Statistics.Jensen_Shannon),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.Statistics.Jensen_Shannon),'$} \\ \hline'];
tobuild{014} = ['\parbox[t]{2mm}{\multirow{7}{*}{\rotatebox[origin=c]{90}{PValues}}} & Kol D & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.PValues.Kolmogorov_D),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.PValues.Kolmogorov_D),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.PValues.Kolmogorov_D),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.PValues.Kolmogorov_D),'$} \\ \cline{2-10}'];
tobuild{015} = ['& CvM & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.PValues.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.PValues.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.PValues.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.PValues.Cramer_von_Mises),'$} \\ \cline{2-10}'];
tobuild{016} = ['& Kuiper & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.PValues.Kuiper),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.PValues.Kuiper),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.PValues.Kuiper),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.PValues.Kuiper),'$} \\ \cline{2-10}'];
tobuild{017} = ['& Watson & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.PValues.Watson),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.PValues.Watson),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.PValues.Watson),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.PValues.Watson),'$} \\ \cline{2-10}'];
tobuild{018} = ['& And-Dar & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.PValues.Anderson_Darling),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.PValues.Anderson_Darling),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.PValues.Anderson_Darling),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.PValues.Anderson_Darling),'$} \\ \cline{2-10}'];
tobuild{019} = ['& Kul-Lei & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.PValues.Kullback_Leibler),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.PValues.Kullback_Leibler),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.PValues.Kullback_Leibler),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.PValues.Kullback_Leibler),'$} \\ \cline{2-10}'];
tobuild{020} = ['& Jen-Sha & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KolD_Pri.PValues.Jensen_Shannon),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.CvM_Pri.PValues.Jensen_Shannon),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Kuiper_Pri.PValues.Jensen_Shannon),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.Watson_Pri.PValues.Jensen_Shannon),'$} \\ \hline\hline'];

tobuild{021} = ['\multicolumn{2}{|c||}{',name,'} & \multicolumn{2}{c||}{\textbf{',labels{5},'}} & \multicolumn{2}{c||}{\textbf{',labels{6},'}} & \multicolumn{2}{c||}{\textbf{',labels{7},'}} & \multicolumn{2}{c||}{\textbf{}} \\ \hline'];
tobuild{022} = ['\multicolumn{2}{|c||}{\textbf{Best Distribution}} & \multicolumn{2}{c||}{\textbf{',distNames{5},'}} & \multicolumn{2}{c||}{\textbf{',distNames{6},'}} & \multicolumn{2}{c||}{\textbf{',distNames{7},'}} & \multicolumn{2}{c||}{\textbf{}} \\ \hline'];
if lowerParaMax == 3
    tobuild{023} = ['\multicolumn{2}{|c||}{\multirow{3}{*}{\textbf{Parameters}}} & ',para.AD_Pri{1},' & ',para.KL_Pri{1},' & ',para.JS_Pri{1},' & ',para.blank{1},'\\ ',clineLower{1}];
    tobuild{024} = ['\multicolumn{2}{|c||}{} & ',para.AD_Pri{2},' & ',para.KL_Pri{2},' & ',para.JS_Pri{2},' & ',para.blank{2},'\\ ',clineLower{2}];
    tobuild{025} = ['\multicolumn{2}{|c||}{} & ',para.AD_Pri{3},' & ',para.KL_Pri{3},' & ',para.JS_Pri{3},' & ',para.blank{3},'\\ ',clineLower{3}];
elseif lowerParaMax == 2
    tobuild{023} = ['\multicolumn{2}{|c||}{\multirow{2}{*}{\textbf{Parameters}}} & ',para.AD_Pri{1},' & ',para.KL_Pri{1},' & ',para.JS_Pri{1},' & ',para.blank{1},'\\ ',clineLower{1}];
    tobuild{024} = ['\multicolumn{2}{|c||}{} & ',para.AD_Pri{2},' & ',para.KL_Pri{2},' & ',para.JS_Pri{2},' & ',para.blank{2},'\\ ',clineLower{2}];
    tobuild{025} = '%This line has been intentionally left blank';    
elseif lowerParaMax == 1
    tobuild{023} = ['\multicolumn{2}{|c||}{\textbf{Parameters}} & ',para.AD_Pri(1),' & ',para.KL_Pri(1),' & ',para.JS_Pri{1},' & ',para.blank{1},'\\ ',clineLower{1}];
    tobuild{024} = '%This line has been intentionally left blank';
    tobuild{025} = '%This line has been intentionally left blank';
end
tobuild{026} = ['\parbox[t]{2mm}{\multirow{7}{*}{\rotatebox[origin=c]{90}{Statistics}}} & Kol D & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.Statistics.Kolmogorov_D),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KL_Pri.Statistics.Kolmogorov_D),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.JS_Pri.Statistics.Kolmogorov_D),'$} & \multicolumn{2}{c||}{} \\ \cline{2-10}'];
tobuild{027} = ['& CvM & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.Statistics.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KL_Pri.Statistics.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.JS_Pri.Statistics.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{} \\ \cline{2-10}'];
tobuild{028} = ['& Kuiper & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.Statistics.Kuiper),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KL_Pri.Statistics.Kuiper),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.JS_Pri.Statistics.Kuiper),'$} & \multicolumn{2}{c||}{} \\ \cline{2-10}'];
tobuild{029} = ['& Watson & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.Statistics.Watson),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KL_Pri.Statistics.Watson),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.JS_Pri.Statistics.Watson),'$} & \multicolumn{2}{c||}{} \\ \cline{2-10}'];
tobuild{030} = ['& And-Dar & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.Statistics.Anderson_Darling),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KL_Pri.Statistics.Anderson_Darling),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.JS_Pri.Statistics.Anderson_Darling),'$} & \multicolumn{2}{c||}{} \\ \cline{2-10}'];
tobuild{031} = ['& Kul-Lei & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.Statistics.Kullback_Leibler),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KL_Pri.Statistics.Kullback_Leibler),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.JS_Pri.Statistics.Kullback_Leibler),'$} & \multicolumn{2}{c||}{} \\ \cline{2-10}'];
tobuild{032} = ['& Jen-Sha & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.Statistics.Jensen_Shannon),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KL_Pri.Statistics.Jensen_Shannon),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.JS_Pri.Statistics.Jensen_Shannon),'$} & \multicolumn{2}{c||}{} \\ \hline'];
tobuild{033} = ['\parbox[t]{2mm}{\multirow{7}{*}{\rotatebox[origin=c]{90}{PValues}}} & Kol D & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.PValues.Kolmogorov_D),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KL_Pri.PValues.Kolmogorov_D),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.JS_Pri.PValues.Kolmogorov_D),'$} & \multicolumn{2}{c||}{} \\ \cline{2-10}'];
tobuild{034} = ['& CvM & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.PValues.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KL_Pri.PValues.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.JS_Pri.PValues.Cramer_von_Mises),'$} & \multicolumn{2}{c||}{} \\ \cline{2-10}'];
tobuild{035} = ['& Kuiper & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.PValues.Kuiper),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KL_Pri.PValues.Kuiper),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.JS_Pri.PValues.Kuiper),'$} & \multicolumn{2}{c||}{} \\ \cline{2-10}'];
tobuild{036} = ['& Watson & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.PValues.Watson),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KL_Pri.PValues.Watson),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.JS_Pri.PValues.Watson),'$} & \multicolumn{2}{c||}{} \\ \cline{2-10}'];
tobuild{037} = ['& And-Dar & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.PValues.Anderson_Darling),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KL_Pri.PValues.Anderson_Darling),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.JS_Pri.PValues.Anderson_Darling),'$} & \multicolumn{2}{c||}{} \\ \cline{2-10}'];
tobuild{038} = ['& Kul-Lei & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.PValues.Kullback_Leibler),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KL_Pri.PValues.Kullback_Leibler),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.JS_Pri.PValues.Kullback_Leibler),'$} & \multicolumn{2}{c||}{} \\ \cline{2-10}'];
tobuild{039} = ['& Jen-Sha & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.AD_Pri.PValues.Jensen_Shannon),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.KL_Pri.PValues.Jensen_Shannon),'$} & \multicolumn{2}{c||}{$',num2matlabstr(SubStruct.JS_Pri.PValues.Jensen_Shannon),'$} & \multicolumn{2}{c||}{} \\ \hline'];


fileID = fopen(filepath,'w');
fprintf(fileID,'%s\r\n',tobuild{:});
fclose(fileID);

end