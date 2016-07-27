function [] = globalAnalysis(inputFolder)

numSplits = 1E2;

timestamp = datestr(now,'yyyymmddTHHMMSS');
iF = ['input/',inputFolder];
toExtract = [iF,'/*.csv'];

oF = ['output_',timestamp];
dir_ref = [oF,'\global'];
mkdir(dir_ref);

fileData = dir(toExtract);
fileList = {fileData.name};

mins = struct('ActiveLinks',Inf(1,6),...
             'InteractionTimes',Inf(1,8),...
             'ActivityPotential',Inf(1,6),...
             'NoContactTimes',Inf(1,6),...
             'NodesActive',Inf(1,6),...
             'Components',Inf(1,6),...
             'Clustering',Inf(1,6),...
             'ComponentNodes',Inf(1,6),...
             'ComponentEdges',Inf(1,6)...
             );
         
maxs = struct('ActiveLinks',-Inf(1,6),...
             'InteractionTimes',-Inf(1,8),...
             'ActivityPotential',-Inf(1,6),...
             'NoContactTimes',-Inf(1,6),...
             'NodesActive',-Inf(1,6),...
             'Components',-Inf(1,6),...
             'Clustering',-Inf(1,6),...
             'ComponentNodes',-Inf(1,6),...
             'ComponentEdges',-Inf(1,6)...
             );         

for i=1:length(fileList)
    currentFile = fileList{i};
    currentClean = strrep(currentFile, '.', '');
    currentClean = strrep(currentClean, '-', '');
    [currentData,current_mins,current_maxs] = analyse(inputFolder,currentFile,'%f %f %f %*s %*s',timestamp);
    %Store Data
    ActiveLinks.(currentClean) = currentData.ActiveLinks_data;
    InteractionTimes.(currentClean) = currentData.InteractionTimes_data;
    ActivityPotential.(currentClean) = currentData.ActivityPotential_data;
    NoContactTimes.(currentClean) = currentData.NoContactTimes_data;
    NodesActive.(currentClean) = currentData.NodesActive_data;
    Components.(currentClean) = currentData.Components_data;
    Clustering.(currentClean) = currentData.Clustering_data;
    ComponentNodes.(currentClean) = currentData.ComponentNodes_data;
    ComponentEdges.(currentClean) = currentData.ComponentEdges_data;
    %Update Mins
    mins.ActiveLinks = min(mins.ActiveLinks,current_mins.ActiveLinks);
    mins.InteractionTimes = min(mins.InteractionTimes,current_mins.InteractionTimes);
    mins.ActivityPotential = min(mins.ActivityPotential,current_mins.ActivityPotential);
    mins.NoContactTimes = min(mins.NoContactTimes,current_mins.NoContactTimes);
    mins.NodesActive = min(mins.NodesActive,current_mins.NodesActive);
    mins.Components = min(mins.Components,current_mins.Components);
    mins.Clustering = min(mins.Clustering,current_mins.Clustering);
    mins.ComponentNodes = min(mins.ComponentNodes,current_mins.ComponentNodes);
    mins.ComponentEdges = min(mins.ComponentEdges,current_mins.ComponentEdges);
    %Update Maxs
    maxs.ActiveLinks = max(maxs.ActiveLinks,current_maxs.ActiveLinks);
    maxs.InteractionTimes = max(maxs.InteractionTimes,current_maxs.InteractionTimes);
    maxs.ActivityPotential = max(maxs.ActivityPotential,current_maxs.ActivityPotential);
    maxs.NoContactTimes = max(maxs.NoContactTimes,current_maxs.NoContactTimes);
    maxs.NodesActive = max(maxs.NodesActive,current_maxs.NodesActive);
    maxs.Components = max(maxs.Components,current_maxs.Components);
    maxs.Clustering = max(maxs.Clustering,current_maxs.Clustering);
    maxs.ComponentNodes = max(maxs.ComponentNodes,current_maxs.ComponentNodes);
    maxs.ComponentEdges = max(maxs.ComponentEdges,current_maxs.ComponentEdges);    
end

Global.ActiveLinks = globalFit_ExpGamRayLN(ActiveLinks,mins.ActiveLinks,maxs.ActiveLinks,numSplits);
Global.InteractionTimes = globalFit_ExpMLGPWei(InteractionTimes,mins.InteractionTimes,maxs.InteractionTimes,numSplits);
Global.ActivityPotential = globalFit_ExpGamRayLN(ActivityPotential,mins.ActivityPotential,maxs.ActivityPotential,numSplits);
Global.NoContactTimes = globalFit_ExpGamRayLN(NoContactTimes,mins.NoContactTimes,maxs.NoContactTimes,numSplits);
Global.NodesActive = globalFit_ExpGamRayLN(NodesActive,mins.NodesActive,maxs.NodesActive,numSplits);
Global.Components = globalFit_ExpGamRayLN(Components,mins.Components,maxs.Components,numSplits);
Global.Clustering = globalFit_ExpGamRayLN(Clustering,mins.Clustering,maxs.Clustering,numSplits);
Global.ComponentNodes = globalFit_ExpGamRayLN(ComponentNodes,mins.ComponentNodes,maxs.ComponentNodes,numSplits);
Global.ComponentEdges = globalFit_ExpGamRayLN(ComponentEdges,mins.ComponentEdges,maxs.ComponentEdges,numSplits);

global_plot(ActiveLinks,Global.ActiveLinks,dir_ref,'ActiveEdges','Fraction of Edges Active',3,1E-3);
global_plot(InteractionTimes,Global.InteractionTimes,dir_ref,'InteractionTimes','Length of Interaction',3,1E-3);
global_plot(ActivityPotential,Global.ActivityPotential,dir_ref,'ActiveEdges','Fraction of Edges Active',3,1E-1);
global_plot(NoContactTimes,Global.NoContactTimes,dir_ref,'TimeBetweenContacts','Length of Time Between Contacts',3,1E-4);
global_plot(NodesActive,Global.NodesActive,dir_ref,'ActiveNodes','Fraction of Nodes Active',1,1E-4);
global_plot(Components,Global.Components,dir_ref,'NumberComponents','Number of Components',0,1E-4);
global_plot(Clustering,Global.Clustering,dir_ref,'GlobalClusteringCoeff','Global Clustering Coefficient',0,1E-2);
global_plot(ComponentNodes,Global.ComponentNodes,dir_ref,'ComponentNodes','Fraction of Nodes per Component',1,1E-1);
global_plot(ComponentEdges,Global.ComponentEdges,dir_ref,'ComponentEdges','Fraction of Edges Active per Component',3,1E-1);

%LATEX