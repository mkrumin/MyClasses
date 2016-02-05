classdef TMazeReplay
    %TMAZEVR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        expRef
        info
        dataTMaze
        timesVRframes
%         timesTrials
%         timesTmazeOnset
%         timesTmazeOffset
        data2p
        times2p
        dataTL
%         timesReward
%         timesSoundStart
        dataBall
        dataEye
        timesEye
        trainingData
        ellData
    end
    
    methods
        %% ========================================================================
        function obj = TMazeReplay(expRef)
            fprintf('Starting the TMazeReplay constructor for experiment %s...\n', expRef);
            obj.expRef = expRef;
            fprintf('Loading and processing ''behavior'' data (which is scrambled Open-Loop replay)...\n');
            obj.dataTMaze = loadTMazeData(obj);
            
            fprintf('Loading 2P data...\n');
            obj.info = ppbox.infoPopulate(expRef);
            for iPlane = 1:obj.info.nPlanes
                fprintf('Plane %d/%d\n', iPlane, obj.info.nPlanes);
                try
                    tmp = load(fullfile(obj.info.folderProcessed, sprintf('%s_plane%03d_ROI', obj.info.basename2p, iPlane)));
                    obj.data2p{iPlane} = tmp.meta;
                    obj.data2p{iPlane}.planeHeaders = obj.data2p{iPlane}.planeHeaders(1);
                catch e
                    fprintf('2-P data for plane %d not loaded\n', iPlane);
                    fprintf(e.message);
                    fprintf('\n');
                end
            end
            clear tmp;
            
            for iPlane = 1:obj.info.nPlanes
                %                 fprintf('Plane %d/%d\n', iPlane, info.nPlanes);
                try
                    if isfield(obj.data2p{iPlane}, 'collageMIP')
                        obj.data2p{iPlane}.avgImage = obj.data2p{iPlane}.collageMIP;
                    else
                        obj.data2p{iPlane}.avgImage = obj.data2p{iPlane}.targetFrame;
                        tiffFileName = fullfile(obj.info.folderProcessed, sprintf('%s_plane%03d_registered_AVG.tiff', obj.info.basename2p, iPlane));
                        if exist(tiffFileName, 'file')
                            tmp = img.loadFrames(tiffFileName, 1, 1, 1);
                            obj.data2p{iPlane}.avgImage = tmp;
                        end
                    end
                catch e
                    fprintf('Average image for plane %d not loaded\n', iPlane);
                    fprintf(e.message);
                    fprintf('\n');
                end
            end
            
            fprintf('Loading Timeline data...\n');
            tmpTL = load(fullfile(obj.info.folderTLLocal, obj.info.basenameTL));
            obj.dataTL = tmpTL.Timeline;
            
            fprintf('Getting ball data...\n');
            obj.dataBall = getRunningSpeed(obj.info);
            
            fprintf('Loading eye-tracking data...\n');
            %             addpath('\\zserver\Code\MouseEyeTrack\');
            try
                eyeFileBasename = dat.expFilePath(obj.expRef, 'eyeTracking');
                if exist([eyeFileBasename{1}, '.mj2'], 'file')
                    eyeFileName = [eyeFileBasename{1}, '.mj2'];
                    obj.dataEye = load([eyeFileBasename{1}, '_processed']);
                else
                    eyeFileName = [eyeFileBasename{2}, '.mj2'];
                    obj.dataEye = load([eyeFileBasename{2}, '_processed']);
                end
                obj.timesEye = et.getFrameTimes(obj.expRef);
            catch
                fprintf('No eye data or eye timestamps data was loaded\n');
            end
            
            fprintf('Getting 2P frame times...\n');
            for iPlane = 1:obj.info.nPlanes
                if ~isempty(obj.data2p{iPlane})
                    fprintf('Plane %d/%d\n', iPlane, obj.info.nPlanes);
                    obj.times2p{iPlane}=ppbox.getFrameTimes(obj.info, obj.data2p{iPlane}.planeFrames);
                end
            end
            
            fprintf('Getting VR stimulus movie times...\n');
            try
                obj.timesVRframes = getVRFrameTimes(obj);
            catch e
                fprintf('Something went wrong...\n');
                fprintf('%s\n', e.message);
            end
            
            % not sure I need these for the scrambled replay analysis,
            % commenting it out for now
            %             fprintf('Getting trials'' times...\n');
            %             obj.timesTrials=ppbox.getStimTimes(obj.info);
            %
            %             nTrials = length(obj.timesVRframes);
            %             obj.timesTmazeOnset = nan(nTrials, 1);
            %             obj.timesTmazeOffset = nan(nTrials, 1);
            %             for iTrial = 2:nTrials
            %                 if ~isempty(obj.timesVRframes(iTrial).t)
            %                     obj.timesTmazeOnset(iTrial) = obj.timesVRframes(iTrial).t(2);
            %                     obj.timesTmazeOffset(iTrial) = obj.timesVRframes(iTrial).t(end);
            %                 end
            %             end
            
            % goodSamples = eyeData.results.goodFit & ~eyeData.results.blink;
            % pupil = interp1(frameTimesEye(goodSamples), eyeData.results.area(goodSamples), ballData.t);
            
            obj.trainingData = cell(obj.info.nPlanes, 1);
            fprintf('\nYay! TMazeReplay constructor fininshed working! \n');
            fprintf('%s dataset is now ready for the analysis!\n\n', obj.expRef);
        end
        %==========================================================================
        function out = Planes(obj)
            out = [];
            for iPlane=1:length(obj.data2p)
                if ~isempty(obj.data2p{iPlane})
                    out = cat(2, out, iPlane);
                end
            end
        end
        
        %==========================================================================
        function [n, attr] = nROIs(obj, iPlane)
            n = size(obj.data2p{iPlane}.F, 2);
            attr = cell2mat(obj.data2p{iPlane}.ROI.CellClasses);
        end
        
        %==========================================================================
        function n = nPixels(obj, iPlane, iROI)
            n = length(obj.data2p{iPlane}.ROI.CellMaps{iROI});
        end

        %==========================================================================
        function n = nTrials(obj)
            n = max([obj.dataTMaze.SESSION.replaySnippets.trialIndex]);
        end

        %==========================================================================
        function ShowROI(obj, iPlane, iROI)
            %             hold off;
            if isfield(obj.data2p{iPlane}, 'avgImage')
                pic = repmat(double(obj.data2p{iPlane}.avgImage), 1, 1, 3);
            elseif isfield(obj.data2p{iPlane}, 'collageMIP')
                obj.data2p{iPlane}.avgImage = obj.data2p{iPlane}.collageMIP;
                pic = repmat(double(obj.data2p{iPlane}.collageMIP), 1, 1, 3);
            else
                obj.data2p{iPlane}.avgImage = obj.data2p{iPlane}.targetFrame;
                pic = repmat(double(obj.data2p{iPlane}.targetFrame), 1, 1, 3);
            end
            
            minMaxValues = prctile(pic(:), [0.1, 99.9]);
            minValue = minMaxValues(1);
            maxValue = minMaxValues(2);
            pic(pic<minValue) = minValue;
            pic(pic>maxValue) = maxValue;
            % the next line is, actually, unnecessary if plotting with imagesc()
            pic = (pic-minValue)/(maxValue-minValue);
            if isequal(size(obj.data2p{iPlane}.avgImage), size(obj.data2p{iPlane}.targetFrame))
                idx = obj.data2p{iPlane}.targetFrameROI{iROI};
            else
                idx = obj.data2p{iPlane}.ROI.CellMaps{iROI};
            end
            [nr, nc, ~] = size(pic);
            pic = reshape(pic, nr*nc, 3);
            pic(idx,1) = 1;
            pic(idx,2) = pic(idx,2)/2;
            pic(idx,3) = pic(idx,3)/2;
            pic = reshape(pic, nr, nc, 3);
            imagesc(pic);
            axis equal tight
            attr = obj.data2p{iPlane}.ROI.CellClasses{iROI};
            expRefString = strrep(obj.expRef, '_', '\_');
            xyz = obj.data2p{iPlane}.ROI.CellXYZMicrons{iROI};
            title({sprintf('%s, Plane #%d, ROI #%d, ''%c''', expRefString, iPlane, iROI, attr); ...
                sprintf('XYZ = (%1.0f,%1.0f,%1.0f) [\\mum], %d pixels', xyz(1), xyz(2), xyz(3), nPixels(obj, iPlane, iROI))});
            axis off;
        end % ShowROI()
        
    end
    
end