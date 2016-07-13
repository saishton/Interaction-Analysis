function [Analysis] = analyse(input_filename,structure,timestamp)

iF = 'input';
oF = ['output_',timestamp];

clean_input = strrep(input_filename, '.', '');
dir_ref = [oF,'\',clean_input];
mkdir(dir_ref);

input = [iF,'/',input_filename];

fid = fopen(input);
rawdata = textscan(fid,structure,'Delimiter',',');
fclose(fid);

%==Extract and Clean Data==%
data = cell2mat(rawdata);
data(:,1) = data(:,1)-data(1,1);
lowestID = min(min(data(:,2)),min(data(:,3)));
data(:,2) = data(:,2)-lowestID+1;
data(:,3) = data(:,3)-lowestID+1;
number_rows = size(data,1);
parfor i=1:number_rows
    thisrow = data(i,:);
    col2 = thisrow(1,2);
    col3 = thisrow(1,3);
    if col2 > col3
        thisrow(1,2) = col3;
        thisrow(1,3) = col2;
        data(i,:) = thisrow;
    end
end
all_IDs = [data(:,2); data(:,3)];
all_active = unique(all_IDs);
num_people = size(all_active,1);
data2 = data(:,2);
data3 = data(:,3);
for i=1:num_people
    oldID = all_active(i);
    data2(data2==oldID) = -i;
    data3(data3==oldID) = -i;
end
data(:,2) = -data2;
data(:,3) = -data3;

%==Perform Analysis==%
[ActiveLinks_FitTool,ActiveLinks_MLE,ActiveLinks_Moments] = analyse_ActiveEdges(data,dir_ref);
[NodesActive_FitTool,NodesActive_MLE,NodesActive_Moments] = analyse_ActiveNodes(data,dir_ref);
[ActivityPotential_FitTool,ActivityPotential_MLE,ActivityPotential_Moments] = analyse_ActivityPotential(data,dir_ref);
[Clustering_FitTool,Clustering_MLE,Clustering_Moments] = analyse_GlobalClusteringCoeff(data,dir_ref);
[InteractionTimes_FitTool,InteractionTimes_MLE,InteractionTimes_Moments] = analyse_InteractionTimes(data,dir_ref);
[Components_FitTool,Components_MLE,Components_Moments] = analyse_NumberComponents(data,dir_ref);
[NoContactTimes_FitTool,NoContactTimes_MLE,NoContactTimes_Moments] = analyse_TimeBetweenContacts(data,dir_ref);
[ComponentNodes_FitTool,ComponentNodes_MLE,ComponentNodes_Moments] = analyse_ComponentNodes(data,dir_ref);
[ComponentEdges_FitTool,ComponentEdges_MLE,ComponentEdges_Moments] = analyse_ComponentEdges(data,dir_ref);

%==Post-Processing & Export==%
datafilename = [dir_ref,'/Distributions.mat'];
save(datafilename,...
    'ActiveLinks_FitTool',...
    'InteractionTimes_FitTool',...
    'ActivityPotential_FitTool',...
    'NoContactTimes_FitTool',...
    'NodesActive_FitTool',...
    'Components_FitTool',...
    'Clustering_FitTool',...
    'ComponentNodes_FitTool',...
    'ComponentEdges_FitTool',...
    'ActiveLinks_Moments',...
    'InteractionTimes_Moments',...
    'ActivityPotential_Moments',...
    'NoContactTimes_Moments',...
    'NodesActive_Moments',...
    'Components_Moments',...
    'Clustering_Moments',...
    'ComponentNodes_Moments',...
    'ComponentEdges_Moments',...
    'ActiveLinks_MLE',...
    'InteractionTimes_MLE',...
    'ActivityPotential_MLE',...
    'NoContactTimes_MLE',...
    'NodesActive_MLE',...
    'Components_MLE',...
    'Clustering_MLE',...
    'ComponentNodes_MLE',...
    'ComponentEdges_MLE'...
)

Analysis = struct(  'ActiveLinks_FitTool',ActiveLinks_FitTool,...
                    'ActiveLinks_Moments',ActiveLinks_Moments,...
                    'ActiveLinks_MLE',ActiveLinks_MLE,...
                    'InteractionTimes_FitTool',InteractionTimes_FitTool,...
                    'InteractionTimes_Moments',InteractionTimes_Moments,...
                    'InteractionTimes_MLE',InteractionTimes_MLE,...
                    'ActivityPotential_FitTool',ActivityPotential_FitTool,...
                    'ActivityPotential_Moments',ActivityPotential_Moments,...
                    'ActivityPotential_MLE',ActivityPotential_MLE,...
                    'NoContactTimes_FitTool',NoContactTimes_FitTool,...
                    'NoContactTimes_Moments',NoContactTimes_Moments,...
                    'NoContactTimes_MLE',NoContactTimes_MLE,...
                    'NodesActive_FitTool',NodesActive_FitTool,...
                    'NodesActive_Moments',NodesActive_Moments,...
                    'NodesActive_MLE',NodesActive_MLE,...
                    'Components_FitTool',Components_FitTool,...
                    'Components_Moments',Components_Moments,...
                    'Components_MLE',Components_MLE,...
                    'Clustering_FitTool',Clustering_FitTool,...
                    'Clustering_Moments',Clustering_Moments,...
                    'Clustering_MLE',Clustering_MLE,...
                    'ComponentNodes_FitTool',ComponentNodes_FitTool,...
                    'ComponentNodes_Moments',ComponentNodes_Moments,...
                    'ComponentNodes_MLE',ComponentNodes_MLE,...
                    'ComponentEdges_FitTool',ComponentEdges_FitTool,...
                    'ComponentEdges_Moments',ComponentEdges_Moments,...
                    'ComponentEdges_MLE',ComponentEdges_MLE,...
                    'NumberPeople',num_people);

analysis2latex(Analysis,dir_ref);
end