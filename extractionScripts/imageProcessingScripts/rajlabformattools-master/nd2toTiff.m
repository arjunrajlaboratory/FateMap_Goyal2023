function [] = nd2toTiff(inputFilesOrFolder,varargin)
% nd2toTiff(inputFilesOrFolder,varargin)
% 
% Reads .nd2 files and produces .tif files
% Creates a .tif for each wavelength and z stack
% Can convert:
%   1 ND2 file --> to a single output directory (1-to-1)
%   many ND2 files --> each to their own directory (N-to-N)
%   many ND2 files --> to a single output directory (N-to-1)
% 
% if inputFilesOrFolder references multiple ND2 files, then a directory
%   will be created for each ND2 file with the same name as the ND2 file by default (N-to-N). However, if
%   a folder name is supplied after the 'outDir' argument, then Tiffs are
%   output to that folder (N-to-1)
% 
% 
% Input
% *|inputFilesOrFolder|*  - filename, a cell array of filenames, or a folder path
% 
% 
% Optional:
%          *|outDir|*  - filepath for the Outputn (Current Dirctory)
%          *|nDigits|* - max order of magnitude for naming convention 
%                           ei. 'dapi0001'  nDigits = 4 (3 )
%          *|channelMap|* - containers.Map object mapping channel names in
%                           nd2 file to names for rajlabimagetools    
%          *|outputXYpositions|* - true or false
% 
% 
% Output
% *|tiff|* - 3D tiff files named "NAME#" were # coresponds to the Z stack and
%           takes into acount previous files. NAME is taken from the wavelength
%           channel (ie DAPI, CY3,...) and changed so that arjlabimagetools can
%           read it.
%           This is an Example of the naming convention as follows:
%                   Our_names     arj_names
%                   alexa594  ->  alexa
%                   Atto647N  ->  cy
%                   cy3       ->  tmr
%                   700       ->  nir
%          Add to this list in the "channelMap"
%
% Required:
%            Get bfmatlab   
%                    1)Go to:   https://downloads.openmicroscopy.org/bio-formats/6.0.1/artifacts/
%                    2)Download bfmatlab.zip
%                    3)Unzip and move bfmatlab folder to your MATLAB folder
%                    4)Add bfmatlab path to Matlab     
%
% Example:
%  >> nd2toTiff('tmr001.nd2');               % read image stack in working directory
% 
% Example:
%  >> nd2toTiff('ABC.nd2'); OR
%  >> nd2toTiff('/Path/to/file/ABC.nd2'); % read from specific filepath, output into the folder 'Path/to/file/ABC/'
%
% Example:
%  >> nd2toTiff('*.nd2');                                       % read all nd2 files from current directory, OR
%  >> nd2toTiff(''/Path/to/folder/');                           % read all nd2 files from specified folder, OR
%  >> inputFilesOrFolder={'FirstFile.nd2';'SecondFile.nd2'};    % specify the nd2 file names
%     nd2toTiff(inputFilesOrFolder);
% 
% Example:
%  >> nd2toTiff(inputFilesOrFolder,"outDir",'/Path/to/output/folder'); % output Tiff files from provided ND2 file(s) into the specified folder
%
% Example:
%  >> nd2toTiff(inputFilesOrFolder,'outputXYpositions',true); % also output a table of the XY coordinates of each position
% 
% Example:
% >> myChannelMap = containers.Map({'Brightfield', 'DAPI', 'YFP', 'YFPa'},...
%                                {'trans'      , 'dapi', 'gfp', 'gfp'});
%    nd2toTiff(inputFilesOrFolder,'channelMap',myChannelMap); % custom  channelMap


% The default MAP
channelMap = containers.Map({'Brightfield', 'DAPI', 'YFP', 'YFP','YFPa','YFPb','YPFc','YPFd', 'CY3','CY3a','CY3b','CY3c','CY3d', 'A594', 'A594a', 'A594b', 'A594c', 'A594d', 'CY5','CY5a','CY5b','CY5c','CY5d', 'A647', '700', 'CY7','NIR'},...
                            {'trans'      , 'dapi', 'gfp', 'gfp','gfpa','gfpb','gfpc','gfpd', 'tmr','tmra','tmrb','tmrc','tmrd', 'alexa','alexaa','alexab','alexac','alexad', 'cy','cya', 'cyb', 'cyc', 'cyd',  'cy'  , 'nir', 'nir','nir'});


% Input check
p = inputParser;
p.addRequired('inputFilesOrFolder')
p.addParameter('outDir', '', @ischar);
p.addParameter('nDigits', 3, @(x)validateattributes(x,{'numeric'},{'positive','integer'}));
p.addParameter('channelMap',channelMap,@(x)isa(channelMap,'containers.Map'));
p.addParameter('outputXYpositions',false,@islogical);

p.parse(inputFilesOrFolder, varargin{:});

channelMap=p.Results.channelMap;
outputXYpositions=p.Results.outputXYpositions;


inputFilesOrFolder = p.Results.inputFilesOrFolder;

%  ------------------------------------------------------------------------
%                   What are the input .nd2 files
%  ------------------------------------------------------------------------
% inputFilesOrFolder = '*.nd2' % current directory
% inputFilesOrFolder = 'path/to/folder'
% inputFilesOrFolder = {'ABC.nd2'; 'XYZ.nd2'};



if ischar(inputFilesOrFolder)
    
    if isfolder(inputFilesOrFolder) % then it's a folder
        if isfile(inputFilesOrFolder), error("input is ambiguous - it is both a folder and a file"),end
        % get nd2 files in this folder
        listOfInputFiles=listFilesWithExt('.nd2','directoryPath',inputFilesOrFolder);
        listOfInputFiles=join([repmat({inputFilesOrFolder},length(listOfInputFiles),1),listOfInputFiles],filesep);
    else
        [a,b,c] = fileparts(inputFilesOrFolder);
        
        if ~strcmp(c,'.nd2'), error(sprintf('input %s must end in .nd2',inputFilesOrFolder)),end
        
        if strcmp(b,'*')
            listOfInputFiles=listFilesWithExt('.nd2');
        elseif isfile(inputFilesOrFolder)
            listOfInputFiles={inputFilesOrFolder};
        else
            error("could not find the input folder or file")
        end
    end
elseif iscell(inputFilesOrFolder)
        
        for i=1:numel(inputFilesOrFolder)
            if ~isfile(inputFilesOrFolder{i})
                error(sprintf('one of the %i filenames provided cannot be found',numel(inputFilesOrFolder)))
            end
        end
        listOfInputFiles=inputFilesOrFolder;
end


%  ------------------------------------------------------------------------
%                   What is the list of output folder(s)
%  ------------------------------------------------------------------------

listOfOutDir=cell('');
if isempty(p.Results.outDir) % default to saving a given nd2 file's tiffs into their own folder (named the same thing as the nd2 file)
    for i=1:numel(listOfInputFiles)
            [a,b,~]=fileparts(listOfInputFiles{i});
            listOfOutDir{i,1} = fullfile(a,b);
    end
else
    listOfOutDir = {p.Results.outDir}; %output to a single folder indicated by the user
    %fprintf('WARNING: All Tiff files will be output to a single folder: %s\n',listOfOutDir{1})
end

fprintf('----------------INPUT-OUTPUT SUMMARY------------------------\n')
fprintf('will read and export these %i nd2 files:\n',numel(listOfInputFiles))
fprintf('   %s\n',listOfInputFiles{:})
fprintf('to the following %i folder(s):\n',numel(listOfOutDir))
fprintf('   %s\n',listOfOutDir{:})
fprintf('------------------------------------------------------------\n\n')

%  ------------------------------------------------------------------------
%                          Read Through nd2  
%  ------------------------------------------------------------------------
nDigits = num2str(p.Results.nDigits);

%cnt_mult_files = imageCount;
for f = 1:numel(listOfInputFiles)
    cnt_stacks = 0; % for this nd2 file
    
    if mod(f-1,1) == 0 %make mod(f-1,5) for every 5th to print
        % file is being read
        fprintf('Working on file %i of %i:\n   %s\n',f,numel(listOfInputFiles),listOfInputFiles{f});
    end
    
    % determine the output directory and whether there are previous .tif
    % files there
    if numel(listOfOutDir)>1
        outDir=listOfOutDir{f};
    else
        outDir=listOfOutDir{1};
    end
    
    if ~isfolder(outDir)
        mkdir(outDir)
    end
    
    % Assign starting number for output files
    if isempty(dir(fullfile(outDir, '*.tif')))
        imageCount = 0;
    else
        currentFiles = dir(fullfile(outDir, '*.tif'));
        currentFiles = {currentFiles.name};
        currentFileCount = regexp(currentFiles, '\d+', 'match');
        currentFileCount = vertcat(currentFileCount{:});
        currentFileCount = str2num(cell2mat(currentFileCount));
        imageCount = max(currentFileCount);
        warning(sprintf('%i tif files currently exist in %s\nlarger filename numbers will be used to avoid overwriting',imageCount,outDir))
    end
    cnt_mult_files = imageCount;
    
    
    % Read through the nd2
    inputFile = listOfInputFiles{f};
    
    reader = bfGetReader(char(inputFile));
    
    omeMeta = reader.getMetadataStore();
    
    Zstacknumb = omeMeta.getImageCount();                    % number of Z stacks
    
    listStackID=nan(Zstacknumb,1);
    listXpositions=nan(Zstacknumb,1);
    listYpositions=nan(Zstacknumb,1);
    
    for i = 1:Zstacknumb % Running through # of z stacks
        cnt_stacks = cnt_stacks + 1;
        
        reader.setSeries(i-1)                                % set ith z stack 
        
        if mod(i,25) == 0
            % stack being read
            fprintf('                       Stack = %03d\n',i);
        end
        
        % Getting info from .nd2 file 
        stackSizeC = omeMeta.getPixelsSizeC(i-1).getValue(); % # of wavelength channels
        stackSizeZ = omeMeta.getPixelsSizeZ(i-1).getValue(); % # of Z slices

        if outputXYpositions
            % Get X and Y positions for later output
            planePositionX=char(omeMeta.getPlanePositionX(i-1,0));
            planePositionX=char(regexp(planePositionX,'\[\S*\]','match'));
            planePositionX=planePositionX(2:end-1);
            planePositionX=str2num(planePositionX);
            
            planePositionY=char(omeMeta.getPlanePositionY(i-1,0));
            planePositionY=char(regexp(planePositionY,'\[\S*\]','match'));
            planePositionY=planePositionY(2:end-1);
            planePositionY=str2num(planePositionY);
            
            
            listXpositions(i)=planePositionX;% look at the first plane in stack. all (from 0 to stackSizeC*stackSizeZ-1) will have same X and Y values
            listYpositions(i)=planePositionY;
            % also save stack ID
            listStackID(i)=cnt_mult_files + cnt_stacks;
            
        end
        % Arrays
        wave_numb = repmat(1:stackSizeC,1,stackSizeZ);     % repeating vector of channel #
        Z_numb    = repmat(1:stackSizeZ,1,stackSizeC);     % repeating vector of Z #

        for ii = 1:stackSizeC  % (# channels in i)
            cnt = 1;
            for iii = 1:stackSizeZ  % (# z slices in stack)
                
                % Read plane from series iSeries at Z, C, T coordinates (iZ, iC, iT)
                iPlane = reader.getIndex(Z_numb(iii) - 1, wave_numb(ii) - 1, 0) + 1; 
                stack_fig  = bfGetPlane(reader, iPlane);

                %---------------------------------------------------------
                %                         Save Tiff 
                %---------------------------------------------------------
                % Change name 
                channelName = omeMeta.getChannelName(i-1, ii-1);
                channelName = channelName.toCharArray';
 
                outBaseName = strcat('%s%0', nDigits, 'd.tif');
                outputFileName = fullfile(outDir, sprintf(outBaseName,channelMap(channelName),cnt_mult_files + cnt_stacks));

                if cnt == 1
                    imwrite(stack_fig, outputFileName)
                    cnt = 1 + cnt;
                else
                    imwrite(stack_fig, outputFileName, 'writemode', 'append')
                    cnt = 1 + cnt;
                end
                
            end
        end            
    end
    
    if outputXYpositions
        % write point metadata to a table
        if numel(listOfOutDir)>1 %N-to-N --> write XYdata now
            outMetaDataTable=table(listStackID,listXpositions,listYpositions,'VariableNames',{'Zstack','Xposition','Yposition'});             
            writetable(outMetaDataTable,fullfile(outDir,'XYdata.csv'))
        else  %only one output directory (N-to-1 or 1-to-1)
            if f==1
                outMetaDataTable=table(listStackID,listXpositions,listYpositions,'VariableNames',{'Zstack','Xposition','Yposition'});            
            else
                outMetaDataTable=[outMetaDataTable;table(listStackID,listXpositions,listYpositions,'VariableNames',{'Zstack','Xposition','Yposition'})];
            end
            
            if f==numel(listOfInputFiles) %time to output the XYdata
                writetable(outMetaDataTable,fullfile(outDir,'XYdata.csv'))
            end
        end
    end
end


end


%% 
function listOfFiles=listFilesWithExt(extension,varargin)
% listFilesWithExt outputs the file names in a directory with a given extension, where listOfFiles is a cell array of character vectors
% 
% Usage
% listOfFiles=listFilesWithExt(extension)
% listOfFiles=listFilesWithExt(extension,'directoryPath',directoryPath)
% listOfFiles=listFilesWithExt(extension,'outputWithExt',false)
% 
% Example:
% 
% listFilesWithExt('.nd2') looks in the current directory and returns
%
%     {'20200922_170759_414__ChannelDAPI,YFP,CY3,A594,CY5_Seq0000.nd2'}
%     {'20200922_222102_814__ChannelDAPI,YFP,CY3,A594,CY5_Seq0000.nd2'}
% 
% listFilesWithExt('.nd2','directoryPath','data/20X','outputWithExt',false) % lists files in directoryPath where the extensions are not included
% 
%     {'20200922_170759_414__ChannelDAPI,YFP,CY3,A594,CY5_Seq0000'}
%     {'20200922_222102_814__ChannelDAPI,YFP,CY3,A594,CY5_Seq0000'}
% 
% listFilesWithExt(extension,'outputWithExt',false)
%


p=inputParser;
p.addParameter('outputWithExt', true, @islogical);
p.addParameter('directoryPath', pwd, @isfolder)
p.parse(varargin{:})

outputWithExt=p.Results.outputWithExt;
directoryPath=p.Results.directoryPath;

if ~isfolder(directoryPath)
    error("could not find directoryPath provided")
end

temp=dir(directoryPath);
dirContents={temp.name}';

assert(ischar(extension))
if ~startsWith(extension,'.')
    extension=['.',extension];
end

listOfFiles=cell('');
number_of_files_with_ext=0;
for ii=1:length(dirContents)
    [~,fileName,fileExt]=fileparts(dirContents{ii});
    if strcmp(fileExt,extension)
        
        number_of_files_with_ext=number_of_files_with_ext+1;
        
    if outputWithExt
        listOfFiles(number_of_files_with_ext,1)={[fileName,extension]};
    else
        listOfFiles(number_of_files_with_ext,1)={fileName};
    end
    end
    
end

end
