function Data=suit_map2surf(input,varargin)
% Data=suit_map2surf(files,varargin);
% Maps a set of volume files onto a the suit surface 
% representation 
% 
% INPUT: 
% Files:        Nx1 Array of filenames that should be mapped 
%               OR 
%               Structure with fields: 
%               V{N}.dim - voxel sizes (x,y,z)
%               V{N}.mat - affine transform matrix
%               V{N}.dat - vol data (x,y,z)
%               OR
%               Nx3 Matrix of foci with each row representing [x,y,z] of
%               foci in atlas space
% OUTPUT:       xxx x N: Matrix of vertex values (xxx) by N,
%               if input were files or structures
%               OR 
%               N x 2: Matrix of coordinates, if input were foci
% VARARGIN: 
%    'space':   SUIT/FSL/SPM: which normalisation routine was used for the
%               images, determine the atlas surface used 
%    'depths':  Vector of Depths at which the line is sampled. 0 refers to the outer surface, 
%               1 to the inner surface. Default are 6 equally spaced positions 
%               along the line: [0 0.2 0.4 0.6 0.8 1]
%    'pial':    File Name of Pial surfac, only necessary if 'space' is
%               undefined
%    'white':   File Name for white surface (if not default) 
%    'stats':   Function that determined how to integrate data over the different sampling depth 
%               'nanmean'     Is default and used for activation data 
%               'mode'        Should be used discrete labels are sampled and the most frequent label is assigned 
%               'minORmax'    Maps the minumum or maximum values, glass brain projection for statistical maps
%    'flat_dir':  Flat file directory


% (c) j.diedrichsen@ucl.ac.uk, 2015
% Edited by Maedbh King (28/02/2017) to include structures as input.
% Advantage of structure is that it removes the need to create a .nii volume
% file to give as input to this function
% -----------------------------------------------------------------
% -----------------------------------------------------------------
%   Set Defaults and deal with user options  
% -----------------------------------------------------------------
global defaults;
if (isempty(defaults))
    spm('Defaults','fmri');
    global defaults;
end; 
space    = ''; 
pial     = 'PIAL_SUIT.coord.gii';
white    = 'WHITE_SUIT.coord.gii';
flat     = 'FLAT.coord.gii'; 
flat_dir = []; 
depths   = [0 0.2 0.4 0.6 0.8 1];

% For volume data
ignore_zeros=1;         % In the volume-data, should zero be set to NaN? Normally yes for beta weights or t-values
stats=@nanmean;         % How are different sampling depth integrated? 

% For foci data 
depth_tolerance = 3;   % How many mm above or below the surface is it allowed to lie? 

vararginoptions(varargin,{'pial','white','space','stats','ignore_zeros','flat_dir','depths'});
numDepths=length(depths);

if ~strmatch('pial', varargin) 
    pial = sprintf('PIAL_%s.coord.gii',space);
end

if ~strmatch('white', varargin)
    white = sprintf('WHITE_%s.coord.gii',space); 
end; 

SCCSid   = '3.1';
SPMid    = spm('FnBanner',mfilename,SCCSid);

% -----------------------------------------------------------------
%   Determine the location of the flat map files
% -----------------------------------------------------------------
if (isempty(flat_dir))
    spmVer=spm('Ver');
    switch(spmVer)
        case {'SPM5','SPM8'}
            suit_dir=fileparts(which('suit_normalize_dartel.m'));
            flat_dir=fullfile(suit_dir,'flatmap');
        case {'SPM12b','SPM12'}
            flat_dir=fullfile(defaults.tbx.dir{1},'suit','flatmap');
        otherwise
            warning(sprintf('Unknown SPM-version: %s',spmVer));
    end;
end;

% -----------------------------------------------------------------
%   Check the atlas space and load the appropriate surface files 
% -----------------------------------------------------------------
try 
    C1=gifti(fullfile(flat_dir,pial));
    C2=gifti(fullfile(flat_dir,white));
    
catch 
    error('unknown space / pial or white file not found');
end; 

c1=C1.vertices;
c2=C2.vertices; 

% -----------------------------------------------------------------
%   If input is image file or structure: Map those. 
% -----------------------------------------------------------------
if (iscell(input) || ischar(input)) || isstruct(input) % MK
% If character array, make into cell
if (ischar(input))
    for i=1:size(input,1);
        V{i}=deblank(input(i,:));
    end;
else 
    V=input; 
end;
A=[];
firstgood=[];
if iscell(V), % check for cell MK
    for i=1:size(V,2);
        try
            A{i}=spm_vol(V{i});
            if (isempty(firstgood))
                firstgood=i;
            end;
        catch
            warning(sprintf('%s could not be opened',V{i}));
            A{i}=[];
        end;
    end;
else % if struct
    A=V; 
end
V=A;

if (isempty(firstgood)) && ~isstruct(input), % MK
    error('none of the images could be opened');
end; 

if (length(ignore_zeros)==1)
    ignore_zeros=ones(length(V),1)*ignore_zeros;
end;

% -----------------------------------------------------------------
%   Loop over the nodes depth and get the linear voxel indices 
% -----------------------------------------------------------------
% For all sampling depth - get the 

if ~isstruct(V)
    dim = V{firstgood}.dim;
    mat = V{firstgood}.mat;
else
    dim = V.dim;
    mat = V.mat;
    V = {V};     %Converting to cell for rest of script
end

for d=1:numDepths % Number of sample points between white and pial surface
    c=(1-depths(d))*c1+depths(d)*c2;
    [i,j,k]=spmj_affine_transform(c(:,1),c(:,2),c(:,3),inv(mat));
    i=round(i); 
    j=round(j); 
    k=round(k); 
    indices(:,d)=i+j*dim(1)+k*dim(1)*dim(2); 
    outside = i<1 | i>dim(1) | j<1 | j>dim(2) | k<1 | k>dim(3); 
    indices(outside,d)=NaN; 
end;

i=find(~isnan(indices));
data=zeros(size(indices))*NaN;

% -----------------------------------------------------------------
%   get the data and apply the statistics over the applicable voxels 
% -----------------------------------------------------------------
for v=1:length(V);
    if (~isempty(V{v}))
        if isfield(V{v},'fname'),% MK 
            X=spm_read_vols(V{v}); 
        else
            X=V{v}.dat; 
        end
        if (ignore_zeros)
            X(X==0)=NaN;
        end;
        data(i)=X(indices(i));
        Data(:,v)=feval(stats,data')';
    else 
        Data(:,v)=nan(size(c1,1),1); 
    end; 
end;

% -----------------------------------------------------------------
%   Map foci onto the flatmap 
% -----------------------------------------------------------------
else
    FL=gifti(fullfile(flat_dir,flat));
    [N,S]=size(input); 
    if ~isnumeric(input) || S~=3 
        error('input needs to be either files names (cell or strarray) or a N x 3 matrix of foci'); 
    end; 
    C=zeros(size(c1,1),3,numDepths); 
    for d=1:numDepths       % Number of sample points between white and pial surface
        C(:,:,d)=(1-depths(d))*c1+depths(d)*c2;
    end; 
    for n=1:N 
        difference = bsxfun(@minus,C,input(n,:));   %
        sqdistance = squeeze(sum(difference.^2,2)); 
        [~,indx]  = min(min(sqdistance,[],2));          % Find the closest line between pial and white 
        Data(n,:) = FL.vertices(indx,:);                         % Take those x,y coordinates 
        vector    = c2(indx,:) - c1(indx,:);  %Determine projection depth  
        x         = input(n,:) - c1(indx,:);        
        Data(n,3) = vector * x' /(vector*vector');      % Normalised projection 
        if Data(n,3)>1 
            dev=(Data(n,3)-1)*norm(vector);             % In terms of mm
        elseif Data(n,3)<0
            dev=-Data(n,3)*norm(vector);             % In terms of mm
        else 
            dev=0; 
        end; 
        if (dev>depth_tolerance) 
            Data(n,:)=NaN; 
        end; 
    end; 
end;
    

