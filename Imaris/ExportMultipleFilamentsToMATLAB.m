%% Instructions
%  Select Filament and send to Matlab as Skel.mat files
%
%  Installation:
%
%  - Copy this file into the XTensions folder in the Imaris installation directory
%  - You will find this function in the Image Processing menu
%
%  <CustomTools>
%      <Menu>
%       <Submenu name="Filament Functions">
%        <Item name="Export Multiple Filaments to MATLAB" icon="Matlab" tooltip="Save Filaments as .mat">
%          <Command>Matlab::ExportMultipleFilamentsToMATLAB(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpFilament">
%          <Item name="Export Mutiple Filaments to MATLAB" icon="Matlab" tooltip="Save Filaments as .mat">
%            <Command>Matlab::ExportMultipleFilamentsToMATLAB(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%
%  Description: Export all filaments contained in selected "Filaments" object as .mat
%  TODO: Convert aRad and aXYZ from real world coordinates back to pixels


%% connect to Imaris Com interface
function ExportMultipleFilamentsToMATLAB(aImarisApplicationID)

if ~isa(aImarisApplicationID, 'COM.Imaris_Application')
    vImarisServer = actxserver('ImarisServer.Server');
    vImarisApplication = vImarisServer.GetObject(aImarisApplicationID);
else
    vImarisApplication = aImarisApplicationID;
end

%% if testing from matlab (comment when running from Imaris)
%vImarisApplication = actxserver('Imaris.Application');
%vImarisApplication.mVisible = true;

%% the user has to create a scene with some spots and surface
vSurpassScene = vImarisApplication.mSurpassScene;
if isequal(vSurpassScene, [])
    msgbox('Please create Surpass scene!')
    return
end

%% get image size and pixel resolution
tmpDataset = vImarisApplication.mDataset; %get the dataset to retrieve size/resolution
xs = tmpDataset.mSizeX; %X size in pixel
ys = tmpDataset.mSizeY; %Y size in pixel
zs = tmpDataset.mSizeZ; %Z size in pixel
xsReal = tmpDataset.mExtendMaxX - tmpDataset.mExtendMinX; %X size in micron
ysReal = tmpDataset.mExtendMaxY - tmpDataset.mExtendMinY; %Y size in micron
zsReal = tmpDataset.mExtendMaxZ - tmpDataset.mExtendMinZ; %Z size in micron
xr = xsReal/xs; %X pixel resolution (usually micron per pixel)
yr = ysReal/ys; %Y pixel resolution (usually micron per pixel)
zr = zsReal/zs; %Z pixel resolution (usually micron per pixel)

%% Make directory of Filaments in surpass scene
cnt = 0;
for vChildIndex = 1:vSurpassScene.GetNumberOfChildren
    if vImarisApplication.mFactory.IsFilament(vSurpassScene.GetChild(vChildIndex - 1))
        cnt = cnt+1;
        vFilaments{cnt} = vSurpassScene.GetChild(vChildIndex - 1);
    end
end

%% Choose which Filaments object to export
vFilamentsCnt = length(vFilaments);
for n= 1:vFilamentsCnt
    vFilamentsName{n} = vFilaments{n}.mName;
end
cellstr = cell2struct(vFilamentsName,{'names'},vFilamentsCnt+2);
str = {cellstr.names};
[vAnswer_yes,~] = listdlg('ListSize',[200 160], 'PromptString','Choose Filament:', 'SelectionMode','single','ListString',str);
aFilament = vFilaments{vAnswer_yes};

aXYZ = aFilament.GetPositionsXYZ;
aRad = aFilament.GetRadii;
aEdges = aFilament.GetEdges;

%% Separate unique filaments contained in aFilament
[C, ia] = setdiff(aEdges(:,1),aEdges(:,2)); % identifies starting point of each filament as non-overlapping sets of filaments
%C contains the numbers (aXYZ orw IDs) of the unique non-overlapping points between filaments, these represents each filament's beginning point (usually cell soma)
%ia contained the row index in the list of edges aEdges where those unique non-overlapping points are. This variable is useful to reconstruct connectivity of the filament

%% Visualize filaments beginning points in imaris as a set of spots
%tmpDiffXYZ = aXYZ(C+1,:); % find coordinates of these unique points
%vSpotsAPosXYZ = tmpDiffXYZ;
%vSpotsARadius = 6*ones(1,numel(tmpDiff));
%vSpotsAPosT = zeros(1,length(tmpDiff));
%vSpotsA = vImarisApplication.mFactory.CreateSpots;
%vSpotsA.Set(vSpotsAPosXYZ, vSpotsAPosT, vSpotsARadius);
%vSpotsA.mName = sprintf('unique');
%vSpotsA.SetColor(0.0, 1.0, 0.0, 0.0);
%vImarisApplication.mSurpassScene.AddChild(vSpotsA);

%% Retrieve connected edges belonging to each filament from the pool aEdges
vProgressDisplay = waitbar(0, 'Splitting each cell into an individual filament');

Skels = [];
for i=1:numel(ia)
    %% Isolate the edges belonging to each individual filament 
    if i ~= numel(ia)
        tmpEdges = aEdges(ia(i):ia(i+1)-1,:);
    else
        tmpEdges = aEdges(ia(i):end,:);
    end

    tmpEdges = tmpEdges +1; % convert indexes in tmpEdges from Imaris (start with 0) to matlab style (start with 1)
    Skels(i).aEdges = tmpEdges;
    Skels(i).aRad = aRad(unique(tmpEdges));   % Extract radius of points only of points in selected edges
    Skels(i).aXYZ = aXYZ(unique(tmpEdges),:); % Extract XYZ coords only of points in selected edges
    Skels(i).aEdges = Skels(i).aEdges - Skels(i).aEdges(1,1) + 1; % shift index edges so that first element starts with 1

    Skels(i).XYZ(:,1) = ceil(Skels(i).aXYZ(:,1)./xr); % Save coordinates in pixel units
    Skels(i).XYZ(:,2) = ceil(Skels(i).aXYZ(:,2)./yr); % Save coordinates in pixel units
    Skels(i).XYZ(:,3) = ceil(Skels(i).aXYZ(:,3)./zr); % Save coordinates in pixel units

    %% March edges to find matching pairs and recreate the branching pattern
    % For each branch first the march from cell body is computer, then
    % all points in common with previously calculated branches are removed
    % this leaves 1 segment per branch without any diplicate
    vNumberOfSpots = length(Skels(i).aRad); % Store total number of points we need to iterate        
    % Find position of terminal and biforcation points in the filament connectivity (aEdges)
    vNumberOfTerminals = 0;
    vTerminals = [];
    vNumberOfForks = 0;
    vForks = [];
    % Start investigating all spots except root (start from position #2 to exclude root)
    for vSpots = 2:vNumberOfSpots
        vEdge = find(Skels(i).aEdges==vSpots);
        % if current edge is a terminal point of the skeleton, it should be listed once
        if length(vEdge)==1
            % disp('found a terminal point');
            vNumberOfTerminals = vNumberOfTerminals + 1;
            vTerminals(vNumberOfTerminals) = vSpots;
        elseif length(vEdge)>2
            % disp('found fork point');
            vNumberOfForks = vNumberOfForks + 1;
            vForks(vNumberOfForks) = vSpots;
        end
    end
    

    % March backwards from each terminal point to the root in orher to find entire branch path
    vPaths=[]; % keeps an ongoing list of points already assigned to a branch
    for vTerminalIndex = 1:vNumberOfTerminals
        vLength = 1;
        vTerminal = vTerminals(vTerminalIndex); % start from terminal point
        vPath = vTerminal;  % add the terminal point to current path
        
        vFound = true;
        while vFound
            % find among edges which one is connect to current terminal
            % vEdge contains the number to wich vTerminal is connected
            % vSide contains which side of the edge is vTerminal in this connection
            [vEdge, vSide] = find(Skels(i).aEdges==vTerminal);
            
            vFound = false;
            for vNeighborIndex = 1:length(vEdge)
                % looks like is marching both directions here
                vNeighbor = Skels(i).aEdges(vEdge(vNeighborIndex), 3-vSide(vNeighborIndex));
                vNeighborSide = 3-vSide(vNeighborIndex); % Store side of the neighbor in the connectivity
                if vNeighbor<vTerminal
                    vNewTerminal = vNeighbor;
                    vFound = true;
                end
            end
            if vFound
                vLength = vLength + 1;
                vTerminal = vNewTerminal;
                vPath(vLength) = vTerminal;
            end
        end
        
        vPath = fliplr(vPath); % Flip path so that is not going terminal->root but root->terminal instead
        vPath = setdiff(vPath, vPaths); % Remove common part with the paths previously calculated
        vPaths = cat(2,vPaths, vPath); % Add current path to paths
        
        Skels(i).branches(vTerminalIndex).XYZ = Skels(i).XYZ(vPath,:);
        Skels(i).branches(vTerminalIndex).Rad = Skels(i).aRad(vPath);
        %Skels(i).branches(vTerminalIndex).Edges = [1:vLength-1;2:vLength]';       
    end


    %% Calculate total length of neurites and store individual segment points into aSeg 
    Skels(i).aSeg = cat(3,Skels(i).aXYZ(Skels(i).aEdges(:,1),:),Skels(i).aXYZ(Skels(i).aEdges(:,2),:));
    Skels(i).Lengths = sqrt((Skels(i).aSeg(:,1,1)-Skels(i).aSeg(:,1,2)).^2 ...
        + (Skels(i).aSeg(:,2,1)-Skels(i).aSeg(:,2,2)).^2 ...
        + (Skels(i).aSeg(:,3,1)-Skels(i).aSeg(:,3,2)).^2);
    Skels(i).TotalLen = sum(Skels(i).Lengths);

    % Calculate soma properties (position, size)
    Skels(i).SomaPtID = 1; % First point is now always the soma
    Skels(i).SomaRad = Skels(i).aRad(Skels(i).SomaPtID);
    Skels(i).SomaArea = pi * (Skels(i).SomaRad)^2;
    
    Skels(i).ResXYZ = [xr yr zr]; % Store image resolution
    Skels(i).ImgSizeXYZ = [xs ys zs]; % Store original image size
    
    % Calculate area of the enclosing polygon (convex hull) on XY projection
    tmpXall = double(Skels(i).aXYZ(:,1));
    tmpYall = double(Skels(i).aXYZ(:,2));
    tmpHull = convhull(tmpXall, tmpYall);
    Skels(i).PolygonAreaNeurites = polygonArea(tmpXall(tmpHull), tmpYall(tmpHull));
    
    % Calculate sholl analysis (step size = 1 micron0
    Skels(i).sholl = shollAnalysis(Skels(i), 1, false);

    % Refresh progress bar
    waitbar(i/numel(ia));
end
close(vProgressDisplay);

%% Get SaveLocation
TPN = uigetdir;
TPN = [TPN filesep];
save([TPN 'Skels.mat'],'Skels');

disp([num2str(numel(Skels)) ' Filaments exported successfully!']);
end