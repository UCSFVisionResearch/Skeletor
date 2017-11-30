%% Instructions
%
%  Send coordinates of Imaris spots to Matlab
%
%  Installation:
%
%  - Copy this file into the XTensions folder in the Imaris installation directory
%  - You will find this function in the Image Processing menu
%
%  <CustomTools>
%      <Menu>
%       <Submenu name="Filament Functions">
%        <Item name="Reflect filament along Z axis" icon="Matlab" tooltip="Reflect filament along Z axis">
%          <Command>Matlab::reflectFilamentAlongZ(%i)</Command>
%        </Item>
%       </Submenu>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpFilament">
%          <Item name="Reflect filament along Z axis" icon="Matlab" tooltip="Reflect filament along Z axis">
%            <Command>Matlab::reflectFilamentAlongZ(%i)</Command>
%          </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
%
%  Description:
%
%   The User chooses which filament to process, coordinates of filament are
%   flipped along the Z axis.
%
%% Connect to Imaris Com interface
function reflectFilamentAlongZ(aImarisApplicationID)

if ~isa(aImarisApplicationID, 'COM.Imaris_Application')
    vImarisServer = actxserver('ImarisServer.Server');
    vImarisApplication = vImarisServer.GetObject(aImarisApplicationID);
else
    vImarisApplication = aImarisApplicationID;
end

%% the user has to create a scene 
vSurpassScene = vImarisApplication.mSurpassScene;
if isequal(vSurpassScene, [])
    msgbox('Please create a Surpass scene!');
    return;
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

%% make directory of Filaments in surpass scene
cnt = 0;
for vChildIndex = 1:vSurpassScene.GetNumberOfChildren
    if vImarisApplication.mFactory.IsFilament(vSurpassScene.GetChild(vChildIndex - 1));
        cnt = cnt+1;
        vFilaments{cnt} = vSurpassScene.GetChild(vChildIndex - 1);
    end
end

%% choose Filament to flip
vFilamentsCnt = length(vFilaments);
for n= 1:vFilamentsCnt;
    vFilamentsName{n} = vFilaments{n}.mName;
end
cellstr = cell2struct(vFilamentsName,{'names'},vFilamentsCnt+2);
str = {cellstr.names};
[vAnswer_yes,ok] = listdlg('ListSize',[200 160], ... 
    'PromptString','Choose Filament:',...
    'SelectionMode','single',...
    'ListString',str);

aFilament = vFilaments{vAnswer_yes};
[aXYZ,aRad,aEdges]=aFilament.Get; % retrieve coordinates of the Filament
aXYZ(:,3) = zsReal - aXYZ(:,3); % invert Z coordinates of the Filament
aFilament.Set(aXYZ,aRad,aEdges); % Push new coordinates back to Imaris

fprintf('Filament Flipped!')
end
