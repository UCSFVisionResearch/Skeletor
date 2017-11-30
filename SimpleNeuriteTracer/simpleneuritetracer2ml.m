%% Convert Fiji's Simple Neurite Tracer skeleton file to matlab
%  created by Luca Della Santina on 2011-10-12
% 
% Skeleton files are GZipped XML .traces, specifications available at:
% http://fiji.lbl.gov/mediawiki/phase3/index.php/Simple_Neurite_Tracer:_.traces_File_Format
%
% Resulting matlab object is a structure:
% traces
% |
% |- calib = (x,y,z) pixel/micron calibration of the skeletonized image
% |- imgsize = (width, height, depth) size of the skeletonized image
% |- totalLength = total dendritic length of the skeleton in micron
% |- branches = individual branches skeletonized, each item is a polygon
%     |- id = unique identifier for the branch
%     |- idParent = id for the parent of this branch (-1 = primary branch)
%     |- idPrimary = id of the primary branch connected to this
%     |- order = branching order (1 = primary branch = originating at soma)


[tmpName tmpPath] = uigetfile('*.traces', 'Select tracing file');
disp('loading skeleton file');

tmpName = gunzip(fullfile(tmpPath, tmpName));

try
    tmpDoc = xmlread(fullfile(tmpPath, char(tmpName))); 
catch
    % if does not work, use the following
    tmpDoc = xmlread(char(tmpName));
end

traces = struct;

tmpTracings = tmpDoc.getDocumentElement;                       % Get the 'tracings' node
tmpEntries = tmpTracings.getChildNodes;

% Load skeleton file
tmpNode = tmpEntries.getFirstChild;
while ~isempty(tmpNode)
    if strcmp(tmpNode.getNodeName, 'samplespacing')
        % Get image calibration
        
        traces.calib.x = str2num(tmpNode.getAttributes.getNamedItem('x').getNodeValue);
        traces.calib.y = str2num(tmpNode.getAttributes.getNamedItem('y').getNodeValue);
        traces.calib.z = str2num(tmpNode.getAttributes.getNamedItem('z').getNodeValue);
    elseif strcmp(tmpNode.getNodeName, 'imagesize')
        % Get image size
        
        traces.imgsize.width = str2num(tmpNode.getAttributes.getNamedItem('width').getNodeValue);
        traces.imgsize.height = str2num(tmpNode.getAttributes.getNamedItem('height').getNodeValue);
        traces.imgsize.depth = str2num(tmpNode.getAttributes.getNamedItem('depth').getNodeValue);
    elseif strcmp(tmpNode.getNodeName, 'path')
        % Get skeleton branches
        
        if tmpNode.hasAttribute('fittedversionof')
            tmpNode = tmpNode.getNextSibling;
            continue; % skip the note if we're dealing with a fitted version
        end
        tmpBranch = struct;

        % General properties of the branch
        tmpBranch.id = str2num(tmpNode.getAttributes.getNamedItem('id').getNodeValue);
        tmpBranch.length = str2num(tmpNode.getAttributes.getNamedItem('reallength').getNodeValue);
        if isempty(tmpNode.getAttributes.getNamedItem('startson'))
            tmpBranch.idParent = -1;
        else
            tmpBranch.idParent = str2num(tmpNode.getAttributes.getNamedItem('startson').getNodeValue);
            % write here code to store branching position from parent dendrite
        end
        
        % Individual points constituting the branch
        tmpBranch.points = ones(1, 7);
        tmpPoints = tmpNode.getChildNodes;
        tmpPoint = tmpPoints.getFirstChild;
        i = 0;
        while ~isempty(tmpPoint)
            if strcmp(tmpPoint.getNodeName, 'point')
                tmpPos = struct;
                tmpPos.x = str2num(tmpPoint.getAttributes.getNamedItem('x').getNodeValue);
                tmpPos.y = str2num(tmpPoint.getAttributes.getNamedItem('y').getNodeValue);
                tmpPos.z = str2num(tmpPoint.getAttributes.getNamedItem('z').getNodeValue);
                tmpPos.xd = str2num(tmpPoint.getAttributes.getNamedItem('xd').getNodeValue);
                tmpPos.yd = str2num(tmpPoint.getAttributes.getNamedItem('yd').getNodeValue);
                tmpPos.zd = str2num(tmpPoint.getAttributes.getNamedItem('zd').getNodeValue);

                tmpPos.r = 0;
                
                i = i+1;
                tmpBranch.points(i, :)= [tmpPos.x, tmpPos.y, tmpPos.z, tmpPos.xd, tmpPos.yd, tmpPos.zd, tmpPos.r];
            end
            tmpPoint = tmpPoint.getNextSibling;
        end

        % Append branch to current cell branches list
        if ~isfield(traces,'branches')
            traces.branches(1) = tmpBranch;
        else
            traces.branches(numel(traces.branches)+1) = tmpBranch;
        end
    end

    tmpNode = tmpNode.getNextSibling;
end

%% Calculate additional parameters for the skeleton

disp('calculating additional skeleton properties');

% Total dendritic length
tmpTotalLen = 0;
for i=1:numel(traces.branches)
    tmpBranch = traces.branches(i);
    tmpTotalLen = tmpTotalLen + tmpBranch.length;
end
traces.totalLength = tmpTotalLen;

% Branching order (1= primary dendrite)
tmpOrd = 1;
tmpOrdIdx = [-1];
tmpNextOrdIdx = [];

while ~isempty(tmpOrdIdx)

    for i=1:numel(traces.branches)
        if ismember(traces.branches(i).idParent, tmpOrdIdx)
            traces.branches(i).order = tmpOrd;                              % Store branching order value
            tmpNextOrdIdx = cat(1, tmpNextOrdIdx, traces.branches(i).id);   % Populate nodes of the next order
        end
    end
    traces.maxOrder = tmpOrd;
    tmpOrd = tmpOrd + 1;
    tmpOrdIdx = tmpNextOrdIdx;
    tmpNextOrdIdx = [];
end

% Primary dendrite generating each branch
tmpPrim = [];
tmpPrimVertex = []; % initial point of primary dendrites (to find soma)

for i=1:numel(traces.branches)
    if traces.branches(i).idParent == -1
        tmpPrim = cat(1,tmpPrim, traces.branches(i).id);
        tmpPrimVertex(numel(tmpPrim), :) = traces.branches(i).points(1,:);
    end
end

for i=1:numel(tmpPrim)
    tmpBranchList = [tmpPrim(i)];
    for j=1:numel(traces.branches)
        if ismember(traces.branches(j).idParent, tmpBranchList) | ...
           ismember(traces.branches(j).id, tmpBranchList)
           
           tmpBranchList = cat(1, tmpBranchList, traces.branches(j).id);
           traces.branches(j).idPrimary = tmpPrim(i);
        end        
    end
end
traces.primaryDendrites = numel(tmpPrim);

% Calculate area of enclosing polyohon (convex hull)
disp('Calculating convex hull area');
tmpXall = []; 
tmpYall = [];
for i=1:numel(traces.branches)
    %plot(traces.branches(i).points(:,1),...
    %     traces.branches(i).points(:,2),...
    %     'Color', [0 0 0]);
    %hold on;
    tmpXall = cat(1, tmpXall, traces.branches(i).points(:,1));
    tmpYall = cat(1, tmpYall, traces.branches(i).points(:,2));
end
tmpHull = convhull(tmpXall, tmpYall);
traces.areaPolygon = polygonArea(tmpXall(tmpHull), tmpYall(tmpHull))*traces.calib.x^2;
%plot(tmpXall(tmpHull), tmpYall(tmpHull), '--','Color', 'magenta');
%axis off;

clear tmp* i j


%% Plot 2D view of the skeleton coloring branches by primary dendrite

tmpH = figure('Name', 'Skeleton 2D projection. Color code = originating primary branch');
set(tmpH, 'Position', [0 traces.imgsize.height/4 traces.imgsize.width*1.25 traces.imgsize.height/2]);
subplot(1,2,1);
hold on;

% identify primary dendrites indexes
tmpPrim = [];
for i=1:numel(traces.branches)
    if traces.branches(i).idParent == -1
        tmpPrim = cat(1,tmpPrim, traces.branches(i).id);
    end
end

% Plot all the branches coloring them by primary dendrite
for i=1:numel(tmpPrim)
    tmpCol = [rand, rand, rand];
    for j=1:numel(traces.branches)
        if traces.branches(j).idPrimary == tmpPrim(i)
           plot(traces.branches(j).points(:,1),...
                traces.imgsize.height - traces.branches(j).points(:,2),...
                'color', tmpCol);
        end        
    end
end

xlim([0 traces.imgsize.height]);
ylim([0 traces.imgsize.width]);
axis off;

subplot(1,2,2);
tmpX = [];
tmpY = [];

for i=1:numel(traces.branches)
    plot(traces.branches(i).points(:,1),...
         traces.branches(i).points(:,2),...
         'Color', 'blue');
     tmpX = cat(1,tmpX, traces.branches(i).points(:,1));
     tmpY = cat(1,tmpY, traces.branches(i).points(:,1));
end

%clear tmp* i j;