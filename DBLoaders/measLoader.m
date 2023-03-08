classdef measLoader <handle
    % Loading set of AO measurements
    % All measurements folders should be in the same project folder
    % All measurements should be named the same with increasing index
    % 
    
    properties
        ao;
        algo;
        vars;
        grid;
        data;

        curData;
    end
    
    methods (Static)
        function uVars = createUserVars()
            uVars.projectPath  = '';
            uVars.measNameTemp = '';
            uVars.numMeas      = '';
            uVars.loadNew = false;

            uVars.lowIdx  = [];
            uVars.highIdx = [];
            
            uVars.phiHighToLow = true;
            uVars.noiseIdxs = [];
        end
    end
    
    methods
        function this = measLoader()
%             this.ao = acoustoOptics();
            this.algo = Algo();
        end
        
        function setVars(this, uVars)
            this.vars.projectPath = uVars.projectPath;
            this.vars.measNameTemp = uVars.measNameTemp;
            this.vars.numMeas = uVars.numMeas;
            this.vars.loadNew = uVars.loadNew;

            this.vars.idxLow = uVars.idxLow;
            this.vars.idxHigh = uVars.idxHigh;
            this.vars.phiHighToLow = uVars.phiHighToLow;

            this.vars.noiseIdxs = uVars.noiseIdxs;
        end
        
        function data = loadMeas(this)
            measNameTemp = sprintf("%s/%s", this.vars.projectPath, this.vars.measNameTemp);
            for i=1:this.vars.numMeas
                filename = sprintf("%s/AO-Results-New.mat", sprintf(measNameTemp, i));
                if ~this.vars.loadNew || ~exist(filename, "file")
                    filename = sprintf("%s/AO-Results.mat", sprintf(measNameTemp, i));
                end
                data.resMeas(i) = load(filename);
                filename = sprintf("%s/AO-Vars.mat", sprintf(measNameTemp, i));
                data.varsMeas(i) = load(filename);
                depthVec = data.varsMeas(i).measVars.algo.len.depthVec*1e3;

                data.phi(i,:)     = data.resMeas(i).phi;
                data.phiRaw(i,:)  = data.resMeas(i).rawPhi;
                data.phiNorm(i,:) = data.resMeas(i).phiNorm;
                data.phiLog(i,:)  = data.resMeas(i).phiLog;
                data.idxPeak(i)   = data.varsMeas(i).measVars.algo.general.extPeakIdx;
                data.depthPeak(i) = depthVec(data.idxPeak(i));
            end
            
            this.data.varsMeas =  data.varsMeas;

            this.data.raw = data;
            this.grid.depthVecRaw = this.data.varsMeas(1).measVars.algo.len.depthVec*1e3;
            this.grid.dDepthRaw   = this.data.varsMeas(1).measVars.algo.len.dDepth*1e3; 

            this.curData.depthVec = this.grid.depthVecRaw;
            this.curData.dDepth   = this.grid.dDepthRaw;
            this.curData.phi      = data.phi; 
            this.curData.phiRaw   = data.phiRaw;
            this.curData.phiNorm  = data.phiNorm;  
            this.curData.phiLog   = data.phiLog; 
            this.curData.idxPeak  = data.idxPeak;
            this.curData.depthPeak  = data.depthPeak;
        end
        
        function resetCurData(this)
            if isfield(this.data, 'cut');      this.data = rmfield(this.data, 'cut'); end
            if isfield(this.data, 'int');      this.data = rmfield(this.data, 'int'); end
            if isfield(this.data, 'align');    this.data = rmfield(this.data, 'align'); end
            if isfield(this.data, 'alignMax'); this.data = rmfield(this.data, 'alignMax'); end

            this.curData.depthVec = this.grid.depthVecRaw;
            this.curData.dDepth   = this.grid.dDepthRaw;
            this.curData.phi      = this.data.raw.phi; 
            this.curData.phiRaw   = this.data.raw.phiRaw;
            this.curData.phiNorm  = this.data.raw.phiNorm;  
            this.curData.phiLog   = this.data.raw.phiLog;
        end

        function cutPhi(this)
            lowIdx  = this.vars.idxLow;
            highIdx = this.vars.idxHigh;
            
            phiCut     = this.curData.phi(:,lowIdx:highIdx);
            phiRawCut  = this.curData.phiRaw(:,lowIdx:highIdx);
            phiNormCut = this.curData.phiNorm(:,lowIdx:highIdx);
            phiLogCut  = this.curData.phiLog(:,lowIdx:highIdx);

            depthVec    = this.data.varsMeas(1).measVars.algo.len.depthZero*1e3;
            depthVecCut = depthVec(lowIdx:highIdx);
            
            this.data.cut.phi     = phiCut;
            this.data.cut.phiRaw  = phiRawCut;
            this.data.cut.phiNorm = phiNormCut;
            this.data.cut.phiLog  = phiLogCut;

            this.grid.depthVecCut = depthVecCut;
            this.grid.dDepthCut   = this.grid.dDepthRaw;

            this.curData.depthVec = this.grid.depthVecCut;
            this.curData.dDepth   = this.grid.dDepthCut;
            this.curData.phi      = phiCut; 
            this.curData.phiRaw   = phiRawCut;
            this.curData.phiNorm  = phiNormCut;  
            this.curData.phiLog   = phiLogCut;
        end
                
        function fixSpeedOfSound(this, c)
            if ~isempty(c)
                this.grid.depthVecFixC = (this.vars.cur.depthVec / this.data.varsMeas(1).measVars.algo.geometry.c) * c;
                this.grid.dDepthFixC   = (this.vars.cur.dDepth / this.data.varsMeas(1).measVars.algo.geometry.c) * c;

                this.curData.depthVec = this.grid.depthVecFixC;
                this.curData.dDepth   = this.grid.dDepthCut;
            end
        end
            
        function interp(this, dDepthInt)
            % Interpolate:
            depthVec    = this.curData.depthVec;
            depthVecInt = depthVec(1): dDepthInt :depthVec(end);
            phiInt      = interp1(depthVec, this.curData.phi',     depthVecInt', 'cubic')';
            phiRawInt   = interp1(depthVec, this.curData.phiRaw',  depthVecInt', 'cubic')';
            phiNormInt  = interp1(depthVec, this.curData.phiNorm', depthVecInt', 'cubic')';
            phiLog = this.curData.phiLog;
            valsLog = mink(unique(phiLog(:)),2);
            phiLog(phiLog == -inf) = valsLog(2);
            phiLogInt  = interp1(depthVec, phiLog',  depthVecInt', 'cubic')';
            
            depthPeak   = this.curData.depthPeak;
            [~, idxPeakInt] = min(abs(depthVecInt' - depthPeak), [], 1);
            depthPeakInt = depthVecInt(idxPeakInt);

            % Collect:
            this.grid.depthVecInt     = depthVecInt;
            this.grid.dDepthInt       = dDepthInt;
            this.data.int.phi         = phiInt;
            this.data.int.phiRaw      = phiRawInt;
            this.data.int.phiNorm     = phiNormInt;
            this.data.int.phiLog      = phiLogInt;
            this.data.int.idxPeak     = idxPeakInt;
            this.data.int.depthPeak   = depthPeakInt;

            this.curData.depthVec = depthVecInt;
            this.curData.dDepth   = dDepthInt;
            this.curData.phi      = phiInt; 
            this.curData.phiRaw   = phiRawInt;
            this.curData.phiNorm  = phiNormInt;  
            this.curData.phiLog   = phiLogInt;
            this.curData.idxPeak  = idxPeakInt;
            this.curData.depthPeak  = depthPeakInt;
        end
        
        function alignDepthVec(this)
            % Align:           
            depthVec = this.curData.depthVec;
            depthVecAligned = depthVec - this.curData.depthPeak(1);

            % Collect:
            this.grid.depthVecAligned = depthVecAligned;
            this.grid.dDepthAligned   = this.curData.dDepth;
            this.data.align.phi       = this.curData.phi;
            this.data.align.phiRaw    = this.curData.phiRaw;
            this.data.align.phiNorm   = this.curData.phiNorm;
            this.data.align.phiLog    = this.curData.phiLog;

            this.curData.depthVec = depthVecAligned;
        end

        function alignToSignalMax(this)
            % Align:
            peakIdx = this.curData.idxPeak;
            offset = peakIdx(1) - peakIdx;

            phiAligned     = zeros(size(this.curData.phi));
            phiRawAligned  = zeros(size(this.curData.phiRaw));
            phiNormAligned = zeros(size(this.curData.phiNorm));
            phiLogAligned  = zeros(size(this.curData.phiLog));

            for i =1:this.vars.numMeas
                phiAligned(i,:)     = circshift(this.curData.phi(i,:), offset(i));
                phiRawAligned(i,:)  = circshift(this.curData.phiRaw(i,:), offset(i));
                phiNormAligned(i,:) = circshift(this.curData.phiNorm(i,:), offset(i));
                phiLogAligned(i,:)  = circshift(this.curData.phiLog(i,:), offset(i));
            end
            
            depthVec = this.curData.depthVec;
            depthVecAligned = depthVec - depthVec(peakIdx(1));

            % Collect:
            this.grid.depthVecAlignSig = depthVecAligned;
            this.grid.dDepthAlignSig   = this.curData.dDepth;
            this.data.alignSig.phi       = phiAligned;
            this.data.alignSig.phiRaw    = phiRawAligned;
            this.data.alignSig.phiNorm   = phiNormAligned;
            this.data.alignSig.phiLog    = phiLogAligned;

            this.curData.depthVec = depthVecAligned;
            this.curData.dDepth   = this.curData.dDepth;
            this.curData.phi      = phiAligned; 
            this.curData.phiRaw   = phiRawAligned;
            this.curData.phiNorm  = phiNormAligned;  
            this.curData.phiLog   = phiLogAligned;
        end

        function normToNoise(this)
            noiseIdxs = this.vars.noiseIdxs;

            % Phi:
            maxVals = max(this.data.phiAligned, [], 2);
            minVals = mean(this.data.phiAligned(:, noiseIdxs),2);
            span = maxVals-minVals;
            this.data.phiAlignNorm = abs((this.data.phiAligned-minVals)./span);
            
            % Phi Raw:
            maxVals = max(this.data.phiRawAligned, [], 2);
            minVals = mean(this.data.phiRawAligned(:,noiseIdxs),2);
            span = maxVals-minVals;
            this.data.phiRawAlignNorm = abs((this.data.phiRawAligned-minVals)./span);
        end

        function newSet = reCalcSet(this, fields)
            for i = 1:this.vars.numMeas
                newSet(i) = this.reCalcMeas(fields, i); 
            end
        end

        function res = reCalcMeas(this, fields, i)
%             this.ao.loadData(this.data.raw.resMeas(i), this.data.raw.varsMeas(i));
            uVars = this.data.raw.varsMeas(i).extVars.algo;
            for j=1:size(fields,1)
                field = fields{j, 1};
                val = fields{j, 2};
                uVars.(field) = val;
            end
            this.algo.setVars(uVars);
            res = this.algo.reconFromFFT(this.data.raw.resMeas(i));
        end

        function resNew = reCalcMeasAndSave(this, fields, i)
            resNew = this.reCalcMeas(fields, i);
            measNameTemp = sprintf("%s/%s", this.vars.projectPath, this.vars.measNameTemp);
            filename = sprintf("%s/AO-Results-New.mat", sprintf(measNameTemp, i));
            save(filename, '-Struct', 'resNew');
        end

        function displayResults(this)
            rows = 1 + isfield(this.data, 'int')+ isfield(this.data, 'align')+ isfield(this.data, 'alignSig');
            cols = 4; 
            gap = 1;
            i=1;
            figure();
            %% Raw Measurement:
            subplot(rows, cols, i)
            plot(this.grid.depthVecRaw, this.data.raw.phi');
            xlabel("Depth [mm]"); ylabel("Fluence [AU]"); title("Phi Meas");
            subplot(rows, cols, i+gap)
            plot(this.grid.depthVecRaw, this.data.raw.phiRaw');
            xlabel("Depth [mm]"); ylabel("Fluence [AU]"); title("Phi Raw");
            subplot(rows, cols, i+2*gap)
            plot(this.grid.depthVecRaw, this.data.raw.phiNorm');
            xlabel("Depth [mm]"); ylabel("Fluence [AU]"); title("Phi Norm");
            subplot(rows, cols, i+3*gap)
            plot(this.grid.depthVecRaw, this.data.raw.phiLog');
            xlabel("Depth [mm]"); ylabel("Fluence [Log]"); title("Phi Log");

            %% Phi Int:
            if isfield(this.data, 'int')
                i = i+cols;
                subplot(rows, cols, i)
                plot(this.grid.depthVecInt, this.data.int.phi')
                xlabel("Depth [mm]"); ylabel("Fluence [AU]"); title("Phi - Interpolated")
                subplot(rows, cols, i+gap)
                plot(this.grid.depthVecInt, this.data.int.phiRaw')
                xlabel("Depth [mm]"); ylabel("Fluence [AU]"); title("Phi Raw - Interpolated")
                subplot(rows, cols, i+2*gap)
                plot(this.grid.depthVecInt, this.data.int.phiNorm')
                xlabel("Depth [mm]"); ylabel("Fluence [AU]"); title("Phi Norm - Interpolated")
                subplot(rows, cols, i+3*gap)
                plot(this.grid.depthVecInt, this.data.int.phiLog')
                xlabel("Depth [mm]"); ylabel("Fluence [Log]"); title("Phi Log - Interpolated")
            end
            
            %% Phi Aligned:
            if isfield(this.data, 'align')
                i = i+cols;
                subplot(rows, cols, i)
                plot(this.grid.depthVecAligned, this.data.align.phi')
                xlabel("Depth [mm]"); ylabel("Fluence [AU]"); title("Phi - Aligned")
                subplot(rows, cols, i+gap)
                plot(this.grid.depthVecAligned, this.data.align.phiRaw')
                xlabel("Depth [mm]"); ylabel("Fluence [AU]"); title("Phi Raw - Aligned")
                subplot(rows, cols, i+2*gap)
                plot(this.grid.depthVecAligned, this.data.align.phiNorm')
                xlabel("Depth [mm]"); ylabel("Fluence [AU]"); title("Phi Norm - Aligned")
                subplot(rows, cols, i+3*gap)
                plot(this.grid.depthVecAligned, this.data.align.phiLog')
                xlabel("Depth [mm]"); ylabel("Fluence [Log]"); title("Phi Log - Aligned")
            end
            
            if isfield(this.data, 'alignSig')
                i = i+cols;
                subplot(rows, cols, i)
                plot(this.grid.depthVecAlignSig, this.data.alignSig.phi')
                xlabel("Depth [mm]"); ylabel("Fluence [AU]"); title("Phi - Aligned")
                subplot(rows, cols, i+gap)
                plot(this.grid.depthVecAlignSig, this.data.alignSig.phiRaw')
                xlabel("Depth [mm]"); ylabel("Fluence [AU]"); title("Phi Raw - Aligned")
                subplot(rows, cols, i+2*gap)
                plot(this.grid.depthVecAlignSig, this.data.alignSig.phiNorm')
                xlabel("Depth [mm]"); ylabel("Fluence [AU]"); title("Phi Norm - Aligned")
                subplot(rows, cols, i+3*gap)
                plot(this.grid.depthVecAlignSig, this.data.alignSig.phiLog')
                xlabel("Depth [mm]"); ylabel("Fluence [Log]"); title("Phi Log - Aligned")
            end
%             %% Phi Cut:
%             if isfield(this.data, 'cut')
%                 i= i+1;
%                 subplot(rows, cols, i);
%                 plot(this.grid.depthVecCut, this.data.cut.phi');
%                 xlabel("Depth [mm]"); ylabel("Fluence [AU]"); title("Phi");
%                 subplot(rows, cols, i+4);
%                 plot(this.grid.depthVecCut, this.data.phiRaw');
%                 xlabel("Depth [mm]"); ylabel("Fluence [AU]"); title("Phi Raw");
%                 subplot(rows, cols, i+8);
%                 plot(this.grid.depthVecCut, this.data.phiNorm');
%                 xlabel("Depth [mm]"); ylabel("Fluence [AU]"); title("Phi Norm");
%                 subplot(rows, cols, i+12);
%                 plot(this.grid.depthVecCut, this.data.phiLog');
%                 xlabel("Depth [mm]"); ylabel("Fluence [Log]"); title("Phi Log");
%             end

%             %% Phi Norm To Noise:
%             if isfield(this.data, 'noiseNorm')
%                 i = i+1;
%                 subplot(rows, cols, i)
%                 plot(this.grid.depthVecIntAligned, this.data.phiAligned')
%                 xlabel("Depth [mm]"); ylabel("Fluence [AU]"); title("Phi - Normalized Noise")
%                 subplot(rows, cols, i+4)
%                 plot(this.grid.depthVecIntAligned, this.data.phiRawAligned')
%                 title("Phi Raw - Normalized Noise")
%                 subplot(rows, cols, i+8)
%                 plot(this.grid.depthVecIntAligned, this.data.phiAligned')
%                 title("Phi Norm - Normalized Noise")
%                 subplot(rows, cols, i+12)
%                 plot(this.grid.depthVecIntAligned, this.data.phiRawAligned')
%                 title("Phi Log- Normalized Noise")
%             end
        end
    end
end

