classdef measLoader <handle
    % Loading set of AO measurements
    % All measurements folders should be in the same project folder
    % All measurements should be named the same with increasing index
    % 
    
    properties
        vars;
        data;
    end
    
    methods (Static)
        function uVars = createUserVars()
            uVars.projectPath  = '';
            uVars.measNameTemp = '';
            uVars.numMeas      = '';
            
            uVars.lowIdx  = [];
            uVars.highIdx = [];
            
            uVars.dX = [];
        end
    end
    
    methods
        function this = measLoader()
            
        end
        
        function setVars(this, uVars)
            this.vars.projectPath = uVars.projectPath;
            this.vars.measNameTemp = uVars.measNameTemp;
            this.vars.numMeas = uVars.numMeas;
            
            this.vars.dX = uVars.dX;
        end
        
        function data = loadMeas(this)
            measNameTemp = sprintf("%s/%s", this.vars.projectPath, this.vars.measNameTemp);
            for i=1:this.vars.numMeas
                filename = sprintf("%s/AO-Results.mat", sprintf(measNameTemp, i));
                data.resMeas(i) = load(filename);
                filename = sprintf("%s/AO-Vars.mat", sprintf(measNameTemp, i));
                data.varsMeas(i) = load(filename);
                
                data.phiMeas(i,:)    = data.resMeas(i).phi;
                data.phiNorm(i,:)    = data.resMeas(i).phiNorm;
                data.phiMeasLog(i,:) = data.resMeas(i).phiLog;
                data.phiRaw(i,:)     = data.resMeas(i).rawPhi;
            end
            
            artfIdx = data.varsMeas(1).measVars.algo.samples.artfctIdxVec;
            data.phiRaw(:,artfIdx) = [];
            this.data = data;
            
            figure();
            subplot(2,2,1)
            plot(data.phiMeas');
            title("Phi Meas");
            subplot(2,2,2)
            plot(data.phiNorm');
            title("Phi Norm");
            subplot(2,2,3)
            plot(data.phiMeasLog');
            title("Phi Log");
            subplot(2,2,4)
            plot(data.phiRaw');
            title("Phi Raw");
        end
        
        function cutPhi(this, lowIdx, highIdx)
            phiCut     = this.data.phiMeas(:,lowIdx:highIdx);
            phiRawCut  = this.data.phiRaw(:,lowIdx:highIdx);
            this.data.phiRaw(:,lowIdx:highIdx);
            
            xVecPreInt = this.data.varsMeas(1).measVars.algo.len.depthZero(1:length(phiCut))*1e3;
            xVecPreIntCut = xVecPreInt(1:size(phiCut,2));
            
            this.data.phiCut         = phiCut;
            this.data.phiRawCut     = phiRawCut;
            
            this.data.xVecRawCut    = xVecPreIntCut;
            
            figure();
            subplot(1,2,1)
            plot(this.data.phiCut')
            subplot(1,2,2)
            plot(this.data.phiRawCut')
        end
                
        function fixSpeedOfSound(this, c)
            this.data.xVecRawCut = (this.data.xVecRawCut / this.data.varsMeas(1).measVars.algo.geometry.c) * c;
        end
            
        function intAndAlign(this, dX)
            xVecInt   = (this.data.xVecRawCut(1)): dX :(this.data.xVecRawCut(end));
            phiInt    = interp1(this.data.xVecRawCut, this.data.phiCut',    xVecInt, 'spline')';
            phiRawInt = interp1(this.data.xVecRawCut, this.data.phiRawCut', xVecInt, 'spline')';
            
            [~, idx] = max(phiInt, [], 2);
            offset = idx(1) - idx;

            phiAligned    = zeros(size(phiInt));
            phiRawAligned = zeros(size(phiRawInt));
            
            for i =1:5
                phiAligned(i,:)    = circshift(phiInt(i,:), offset(i));
                phiRawAligned(i,:) = circshift(phiRawInt(i,:), offset(i));
            end
            
            this.data.xVecInt   = xVecInt;
            this.data.phiInt    = phiInt;
            this.data.phiRawInt = phiRawInt;
            
            this.data.xVecAligned   = xVecInt - xVecInt(idx(1));
            this.data.phiAligned    = phiAligned;
            this.data.phiRawAligned = phiRawAligned;
            
            figure();
            subplot(2,2,1)
            plot(this.data.xVecAligned, this.data.phiAligned')
            subplot(2,2,2)
            plot(this.data.xVecAligned, this.data.phiRawAligned')
            subplot(2,2,3)
            plot( this.data.phiAligned')
            subplot(2,2,4)
            plot(this.data.phiRawAligned')
        end
        
        function normToTail(this, tailIdx)
            % Phi:
            maxVals = max(this.data.phiAligned, [], 2);
            minVals = mean(this.data.phiAligned(:,tailIdx:end),2);
            span = maxVals-minVals;
            this.data.phiAlignNorm = abs((this.data.phiAligned-minVals)./span);
            
            % Phi Raw:
            maxVals = max(this.data.phiRawAligned, [], 2);
            minVals = mean(this.data.phiRawAligned(:,tailIdx:end),2);
            span = maxVals-minVals;
            this.data.phiRawAlignNorm = abs((this.data.phiRawAligned-minVals)./span);

            figure();
            subplot(1,2,1)
            plot(this.data.xVecAligned, this.data.phiAlignNorm')
            subplot(1,2,2)
            plot(this.data.xVecAligned, this.data.phiRawAlignNorm')
            
        end

    end
end

