% E158rep2_tissue_raw data processing
%% Load Tcond
%clear
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/paper';
Tcond=readtable(fullfile(parentDir,'experiments/E158rep2_tissue_Tcond.xlsx'));
%Tcond.Nd2_dir=replace(Tcond.Nd2_dir,"/Volumes/IAND_04/","/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/rawData/");

opts = detectImportOptions(fullfile(parentDir,'experiments/E158rep2_tissue_Tscan.xlsx')); opts = setvartype(opts,{'condIDs'},{'char'});
Tscan=readtable(fullfile(parentDir,'experiments/E158rep2_tissue_Tscan.xlsx'),opts);
Tscan.condIDs=cellfun(@(x) str2num(x),Tscan.condIDs,'UniformOutput',false);
cd(parentDir)

scanIDsList=Tscan.scanID;
%scanIDsList=[3;4]
%scanIDsList=(8)'
%scanIDsList=(6)'
%warning('scan ID list is shorter')

%scanIDsList=1;
%subscanGenerationMode='CornerSubscansOnly';

%% make ScanInfo1.mat file with all rounds, including the XY positions of each tile
% % still need to fix stitching issue - did I shift images in wrong direction.
for iScanID=1:length(scanIDsList)
    
    scanID=scanIDsList(iScanID);
    rowInTscan=find(Tscan.scanID==scanID);
    assert(numel(rowInTscan)==1)
    
    condIDsAll=Tscan.condIDs{rowInTscan};
    numRounds=length(condIDsAll);
    Nd2FilepathList=cell(numRounds,1);
    %channelLabels=repmat({{''}},numRounds,1);
    Tchan=table();
    for thisRound=1:length(Nd2FilepathList)
        
        % make Nd2FilepathList
        rowInTcond=find(Tcond.condID==condIDsAll(thisRound));
        [Nd2_filepath,Nd2_filename]=dirAndFilenametext2filepath(Tcond.Nd2_dir{rowInTcond},Tcond.Nd2_filenametext{rowInTcond});
        Nd2FilepathList{thisRound}=Nd2_filepath;
         % make channelLabels input
        channelLabels=Tcond{rowInTcond, {'A488_gene',   'A555_gene',   'A594_gene', 'A647_gene'}};
        channels=                       {'YFP',         'CY3',         'A594',      'CY5'};
        
        %=~ismember(channelLabels,{''}); channels=channels(idxNotEmpty); channelLabels=channelLabels(idxNotEmpty);
        channels=       [{'DAPI','Brightfield'},channels]; channelLabels=  [{'',''},channelLabels]; % append DAPI and Brightfield
        channelPrefixes=join([repmat(Tcond.cyclePrefix(rowInTcond)   ,1,length(channels));channels],'_',1);
        Tchan=[Tchan; [array2table(repmat(thisRound,length(channels),1),'VariableNames',{'round'}), array2table([channels',channelLabels',channelPrefixes'],'VariableNames',{'channel','channelLabels','channelPrefixes'})]]; % add to channel labels table
    end
    %Tchan
    
    % make R4 (with YFP for cyto segmentations) the reference round
    refRound=1;
    
    %% call makeScanObject
    
    Scan=makeScanObject(Nd2FilepathList,...
       'channelLabels',Tchan,...
       'referenceRound',refRound,...
       'numControlPointsPerRegistration',15);%,...
       %'rowDimIsFlippedVsYDim',true,...
       %'cameraAngleOfReferenceRound',-179.7526);%from wellA1, refRound=1, 8pairs: -179.7526 although was using -179.792 for a while % take out 'cameraAngle' pre-specification later on
    %warning('CAMERA ANGLE FIXED')
    
    Scan.scanName=Tscan.scanName{rowInTscan};
    
    scanDir=Tscan.scanDir{rowInTscan};
    if ~isfolder(scanDir)
        mkdir(scanDir)
    end
    save(fullfile(scanDir,filesep,'ScanObject.mat'),'Scan', '-v7.3')
    
    %     % save Scan
    %     scanDir=Tscan.scanDir{rowInTscan};
    %     if ~isdir(scanDir)
    %         mkdir(scanDir)
    %     end
    %     save(fullfile(scanDir,filesep,'ScanInfo1'),'Scan', '-v7.3')
    
end % end scans loop

%% makeStitches
%scanIDsList=Tscan.scanID
scanIDsList=[2;3;4]
%StitchSubregionSubsetList=[1 1;3 10];
%warning('subregion subset only');
%[a,b]=ind2sub([10,10],1:100);
%SubregionsAll=[b',a'];

%StitchSubregionSubsetList=SubregionsAll(1:16,:);
%StitchSubregionSubsetList=SubregionsAll(17:end,:);

%clear Scan
for iScanID=1:length(scanIDsList)
    scanID=scanIDsList(iScanID);
    fprintf('%s\n',repmat('-',1,200))
    fprintf('%s starting makeStitches for scanID=%i at %s  %s\n',repmat('-',1,80),scanID,char(datetime),repmat('-',1,80))
    
    rowInTscan=find(Tscan.scanID==scanID);
    assert(numel(rowInTscan)==1)

	%SubregionArray=[6 6]; % well B3, with resistant cells, had smaller scan (24x24 instead of 39x39)
    %SubregionArray=[10 10];
    SubregionArray=[Tscan.SubregionRows(rowInTscan),Tscan.SubregionColumns(rowInTscan)];

    
    
    scanDir=Tscan.scanDir{rowInTscan};
    load(fullfile(scanDir,filesep,'ScanObject.mat'))
    
    tifNameFormat={'channelPrefixes','_','channelExposureTimesMs','ms_','channelLabels'};
%     if ismember(scanID,[7,8])
%         tifNameFormat={'channelPrefixes','_','channelExposureTimesMs','ms_','channelLabels'};
%     else
%         tifNameFormat={'R','round','_','channelNames','_','channelLabels'};
%     end
    
    ScanBounds=Tscan.ScanBounds{rowInTscan};
    
    if scanID==2
    makeStitches(Scan,...
        'SubregionArray',SubregionArray,...
        'tifNameFormat',tifNameFormat,...
        'tiffOutParentDir',scanDir,...
        'ScanBounds',ScanBounds,...
        'StitchSubregionSubsetOnly',[1,4])
    else
            makeStitches(Scan,...
        'SubregionArray',SubregionArray,...
        'tifNameFormat',tifNameFormat,...
        'tiffOutParentDir',scanDir,...
        'ScanBounds',ScanBounds)
    end
    %'StitchSubregionSubsetOnly',[2,2])
    
end

