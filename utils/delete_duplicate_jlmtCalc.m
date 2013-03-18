% Maintenance script to delete older ipem calculations with specified duplicate field values
% Analysis files that have illegal data formats (not ensemble data
% struct, missing params data field) or are missing
% required fields will also be deleted.
%
% INSTRUCTIONS
% 1) Set REPORT_ONLY to 0 if you just want a report and no
%    deletions. If you want to delete analyses that meet deletion
%    criteria, then set this to 1.
% 2) Set the log file to a file path and name where you want the
%    report to be logged.
% 3) Set rootDir to the root directory under which the analyses
%    reside.
% 4) Set the field values of fieldMatch to cell arrays of
%    fieldnames which you would like to match. fieldMatch.ani, for
%    example, is a cell array of fieldnames to match for ani
%    calculations. All calculations for which these param value match
%    will be considered duplicate calculations. Only the newest
%    of the duplicate calculations will be saved. Also,
%    calculations that are missing any of these fields will be
%    regarded "illegal" or outdated and will also be deleted.
% 5) Set processCalc to a cell array of calculation names that you
%    want to process. Those not listed here will be skipped. For
%    example, processCalc={'ani'} will only process ani calculations
%    and other calculations such as periodicity pitch or context image
%    calculations will not be checked.
% 6) Run the script. This should take a little while.
%
% Copyright (c) 2008-2013 The Regents of the University of California
% All Rights Reserved
%
% Author(s):
% 10/7/08 Stefan Tomic - First Version

%specify whether you only want a report. If false, deletions will
%occur in addition to report.
REPORT_ONLY = 0;

%file where report will be logged
logFileName = '/tmp/duplicate_ipem_calcs2.txt';

%root directory where analysis mat files will be searched
rootDir = '/data2/stimuli';

%note that you can only specify top level parameters in these
%lists. Lower level parameters will still be compared, by
%specifying the top level field under which they belong. (e.g. if
%you specify 'resonatorBanks' for rp, then rp.resonatorBank.method,
%rp.resonatorBank.spacingFactor, etc. will be included).
fieldMatch.ani = {'normalize_wav','DownSampleFactor','NumOfChannels','FirstCBU','CBUStep','normalize_maxVal'};
fieldMatch.pp = {'LowFreqency','FrameWidth','FrameStepSize','Atten'};
fieldMatch.ci = {'PeriodicityPitch','SnapShot','HalfDecayChords','HalfDecayToneCenters','Enlargement'};
fieldMatch.toract = {'ci','calc_sphere_harm','spher_harm','norm','som','SampleFreq'};
fieldMatch.rp  = {'calcDiff','sumBandsAfter', ... %'mulawCompression'
		  'resonatorBank','inDataType','sumAdjacentBands','ipemRMSMethod',...
		  'Fs','rmsFrameWidth','bypassOnsetDetect','convertInputToMono'};

%calculations to process and check for consistency. Those not specified here will be skipped
%processCalc = {'ani','pp','ci','toract','rp'};
processCalc = {'rp'};

%first column is the directory name in which the calculations
%reside. Second column is the variable name in which they are stored.
calcInfo = {'ani','ani';...
	     'pp','pp';...
	     'ci','cntxt_img';...
	     'toract','toract';...
	     'rp','rp'};

%picking up default params for all the calculations
%this is only used to determine what the field names of the
%substructures (e.g. rp.resonatorBank.method,rp.resonatorBank.spacingFactor) are.
def.ani = calc_ani('getDefaultParams');
def.pp = params_pp;
def.ci = params_ci;
def.toract.ci.siglist = {'sig1','sig2'};
def.toract.calc_spher_harm = 1;
def.toract.spher_harm.nharm_theta = 3;
def.toract.spher_harm.nharm_phi = 4;
def.toract.spher_harm.min_rsqr = 0.95;
def.toract.norm = [];
def.toract.som.fname = [];
def.toract.SampleFreq = [];
def.rp = rhythm_profiler('getDefaultParams');

try
  logFileID = fopen(logFileName,'w');
catch
  disp('Can''t open log file or missing log file name. Outputting messages to stdout.\n\n');
  logFileID = 1;
end

fprintf(['\nObtaining complete list of analysis mat files.\nThis may' ...
	 ' take several minutes...\n']);
allAnFiles = dirrec(rootDir,'.mat');
fprintf('Found %d total analysis files, now parsing through them...\n\n',length(allAnFiles));
fprintf(['Calculations of type: %s will be searched and inspected '...
        'for integrity and parameter matching.\n\n'],cell2str(processCalc,', '));

if(~REPORT_ONLY)
  fprintf(['DELETION FLAG IS SET. Analyses with duplicate parameters' ...
	   ' or illegal data struct/params struct will be deleted.\n\n']);
else
  fprintf(['Deletion flag IS NOT set. Only a report will be generated.' ...
	   ' No files will be deleted.\n\n']);

end

calcTypes = calcInfo(:,1);
nTypes = length(calcTypes);
for iType = 1:length(processCalc)
  thisType = processCalc{iType};
  calcIdx = strmatch(thisType,calcInfo(:,1));  
  varName = calcInfo{calcIdx,2};

  fprintf(logFileID,'FINDING %s ANALYSIS FILES...\n\n',thisType); 
  dirNames = regexp(allAnFiles,['(?<dirName>.*/' thisType ')/.*.mat'],'names');
  %convert dirNames from a cell array to a struct array
  dirNames = [dirNames{:}];
  %now convert to a cell array of strings
  dirNames = {dirNames(:).dirName};
  %find unique list of analysis directories
  dirNames = unique(dirNames);
  
  nDirs = length(dirNames);
  fprintf('Found %d directories with %s analyses.\n',nDirs,thisType);  
  for iDir = 1:nDirs

    fprintf('Processing ''%s'' directory %d of %d...\n',thisType,iDir,nDirs); 
    thisDir = dirNames{iDir};
    fileIdxs = strmatch(thisDir,allAnFiles);
    filesThisDir = allAnFiles(fileIdxs);
    nFiles = length(filesThisDir);
    matchVector = NaN(1,nFiles);
    matchNum = 1;
    for iFile = 1:nFiles
      if(~isnan(matchVector(iFile)))
	%if matchVector for this file is not nan, then we have already matched
        %this file to another one.
	continue
      end
      loadAn = load(filesThisDir{iFile});
      anStruct = loadAn.(varName);
      %if this is not a valid ensemble data struct, or there is no
      %'params' data, then the 'params' cell for this file will be
      %empty. We can choose what to do with these later.
      if(isfield(anStruct,'data') && isfield(anStruct,'vars'))
	anStructCols = set_var_col_const(anStruct.vars);
	if(isfield(anStructCols,'params'))
	  fileParams = anStruct.data{anStructCols.params};

	  %parameters may be a structure contained in a single cell. This
          %may (should) be an obsolete format, but the data may be
          %good, so we're going to check these normally
	  if(iscell(fileParams))
	    fileParams = fileParams{1};
	  end
	  
	  pTypes = fieldnames(fileParams);
	  %go through the params for each calc type contained in a
          %particular file ('rp' calc params will also have 'ani' calc params, etc.)
	  for ipType = 1:length(pTypes)
	    thisPType = pTypes{ipType};
	    
	    %determine if the file has all the parameters (at top level
            %as well as substructures). If not, they are incomplete
            %and probably outdated
	    defSubStruct = def.(thisPType);
	    removeFieldsFromDef = setdiff(fieldnames(defSubStruct),fieldMatch.(thisPType));
	    defSubStruct = rmfield(defSubStruct,removeFieldsFromDef);
	    %perform the match
	    subStructMatch = compare_structs(defSubStruct,fileParams.(thisPType),'values',0,'substruct',1,'types',0);
	    if(~subStructMatch)
	      matchVector(iFile) = -1;
	      break; %no need for further checking, break out of inner loop
	    end
	    
	    %remove fields from our file parameters that we won't
            %be comparing 
	    paramFieldNames = fieldnames(fileParams.(thisPType));
	    removeFieldNames = setdiff(paramFieldNames,fieldMatch.(thisPType));
	    fileParams.(thisPType) = rmfield(fileParams.(thisPType),removeFieldNames);
	  end

	  %if we have determined the status of this file at this
          %point (e.g. incomplete params), move to the next file
	  if(~isnan(matchVector(iFile)))
	    continue;
	  end
	  
	  %find all files with the same params
	  %we are only checking subsets of params (those that were
          %specified in fileParams)
	  cmpParams.subset = 1;
	  cmpParams.allMatches = 1;
	  cmpParams.paramFind = fileParams;
	  cmpParams.varName = varName;
	  matchingFiles = check_anal_exist(thisDir,cmpParams);
	  [pth,thisFstub,ext] = fileparts(filesThisDir{iFile});
	  matchingFiles = setdiff(matchingFiles,[thisFstub ext]);

	  matchVector(iFile) = matchNum;
	  if(~isempty(matchingFiles))
	    matchingFpaths = strcat([thisDir '/'],matchingFiles);
	    [foundMatches,matchingIdxs] = ismember(matchingFpaths,filesThisDir);
	    matchVector(matchingIdxs) = matchNum;
	  end
	  matchNum = matchNum + 1;
	else
	  matchVector(iFile) = 0;	  
	end %if(isfield...
      else
	matchVector(iFile) = 0;
      end %if(isfield...
    end%for iFile...
  
    idxIllegal = find(matchVector == 0);
    nIllegal = length(idxIllegal);
    if(nIllegal > 0)
      fprintf(logFileID,['Found %d file(s) with illegal or missing' ...
			 ' data struct or absent params struct:\n'],nIllegal);
      fprintf(logFileID,'%s\n',filesThisDir{idxIllegal});
      
      if(~REPORT_ONLY)
	for iIll = idxIllegal
	  fprintf(logFileID,'Deleting these files...\n\n');
	  delete(filesThisDir{iIll});
	end
      else
	fprintf(logFileID,'\n\n');
      end
    end
      
    idxMissing = find(matchVector == -1);
    nMissing = length(idxMissing);
    if(nMissing > 0)
      fprintf(logFileID,['Found %d file(s) that are missing individual parameters.' ...
			 ' These files are probably outdated.\n'],nMissing);
      fprintf(logFileID,'%s\n',filesThisDir{idxMissing});
      
      if(~REPORT_ONLY)
	for iM = idxMissing
	  fprintf(logFileID,'Deleting these files...\n\n');
	  delete(filesThisDir{iM});
	end
      else
	fprintf(logFileID,'\n\n');
      end
    end

    
    %a group is a set of files which have a match on all their parameters
    nGroups = max(matchVector);
    if(nGroups < 1)
      continue
    end
    
    %find the newest file in each group, either perform a report or deletion
    dirInfo = dir(thisDir);
    dirInfoFiles = {dirInfo.name};
    dirInfoDateStr = {dirInfo.date};
    dirInfoDateNum = datenum(dirInfoDateStr,0);
    
 
 

    for iGroup = 1:nGroups
      grpIdxs = find(matchVector == iGroup);
      grpFiles = filesThisDir(grpIdxs);
      [dummy,fileMatchIdx] = ismember(grpFiles,strcat(thisDir,'/',dirInfoFiles));
      grpDates = dirInfoDateNum(fileMatchIdx);
      [dummy,newestGrpIdx] = max(grpDates);
      
      if(length(grpIdxs) == 1)
	fprintf(logFileID,'Found %d unique file(s): \n',length(grpIdxs));
      else
	fprintf(logFileID,'Found %d similar file(s): \n',length(grpIdxs));
      end
	
      reportFnames = grpFiles;
      reportFnames{newestGrpIdx} = [reportFnames{newestGrpIdx} ...
		    ' (newest)'];
      fprintf(logFileID,'%s\n',reportFnames{:});
          
      if(~REPORT_ONLY && (length(grpIdxs) > 1) )
	fprintf(logFileID,'Deleting older files....\n\n');
	delIdxs = grpIdxs;
	delIdxs(newestGrpIdx) = [];
	for dIdx = delIdxs
	  delete(filesThisDir{dIdx});
	end
      else
	fprintf(logFileID,'\n\n');
      end
      
    end
    
  end%for iDir...
end%for iType...


if(logFileID > 1)
  fclose(logFileID);
end
