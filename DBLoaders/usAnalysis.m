classdef usAnalysis < handle
    % US analysis - receive a data file of US transducer pressure field and
    % returns its parameters: 
    % (1) AC coupled Data
    % (2) P2P 3D map
    % (3) speed of sound
    % (4) calculates time delay and fixes the US axis and the time axis
    % (5) focal length
    % (6) 3dB width on transversal and longitudinal axes.
    %
    % This code assumes that:
    % cont axis  = Y
    % disc1 axis = Z
    % disc2 axis = X
    % US axis    = X
    % US near field as at low X
    
    properties
        tf
        vars
        us
        usVars
    end
    
    methods (Static)
        function uVars = createUserVars()
            uVars.usDataPath = [];
            uVars.usDataType = [];
            uVars.intFactor  = [];
        end
    end
    
    methods
        function this = usAnalysis()
            this.usVars.curDataPath = '';
        end
        
        function setVars(this, uVars)
           this.vars.usDataPath = uVars.usDataPath;
           this.vars.usDataType = uVars.usDataType; 
           this.vars.intFactor = uVars.intFactor;
        end
        
        function vars = getVars(this)
            vars = this.usVars;
        end

        function [data, vars] = analyse(this)
            if ~strcmp(this.usVars.curDataPath, this.vars.usDataPath)
                
                this.loadUS();
                this.acCoupling();
                this.extractP2P();
                this.convertToPascal();
                this.interpolate();
                this.extractFocusIdx();
                this.extractProfiles();
                this.extractSpeedOfSound();
                this.calcDepthVec();
                this.calcMinimalDelay()
                this.fixUSAxisOffset();
                this.calcFocalLen();
                this.calcRayleighLen();
                this.calcWaistSize();
                this.calcFocalFFT();
                this.calcEnvelopes();
                this.cutFocalData();
                this.interpolateFocalEnvelopes();
                this.createSpatialData();
                this.calcPulseWidth();
                this.cutPulses();
                this.calaPulsesProfile();
                this.displayResults();
            end

            data = this.us;
            vars = this.usVars;
            
        end
        
        function loadUS(this)
            fprintf("US-Analysis: Loading \n");
            if ~strcmp(this.usVars.curDataPath, this.vars.usDataPath)
                this.us = [];
                res = load(this.vars.usDataPath);
                res.resCs = flip(res.resCs, 3); % align X
                this.us.data = squeeze(res.resCs); %remove channel
                this.usVars.raw = res.csVars;
                this.usVars.rawGenVars = res.genVars;
                this.usVars.curDataPath = this.vars.usDataPath;
                this.usVars.calcEnv = true;
            end
            
            switch this.vars.usDataType
                case '3D'
                    this.usVars.dxRaw = this.usVars.raw.axDisc2Stride; %[mm]
                    this.usVars.dyRaw = this.usVars.raw.axScanStride;  %[mm]
                    this.usVars.dzRaw = this.usVars.raw.axDisc1Stride; %[mm]

                    this.usVars.xVecRaw = this.usVars.raw.disc2Vec; %[mm]
                    this.usVars.yVecRaw = this.usVars.raw.scanVecBin - mean(this.usVars.raw.scanVecBin);  %[mm]
                    this.usVars.zVecRaw = this.usVars.raw.disc1Vec- mean(this.usVars.raw.disc1Vec); %[mm]

                    this.usVars.tVec = this.usVars.raw.tVec;
                case '2D'
                    this.usVars.dxRaw = this.usVars.raw.axDisc2Stride; %[mm]
                    this.usVars.dyRaw = this.usVars.raw.axScanStride;  %[mm]
                    this.usVars.dzRaw = this.usVars.raw.axDisc1Stride; %[mm]

                    this.usVars.xVecRaw = this.usVars.raw.disc2Vec; %[mm]
                    this.usVars.yVecRaw = this.usVars.raw.scanVecBin - mean(this.usVars.raw.scanVecBin);  %[mm]
                    this.usVars.zVecRaw = this.usVars.raw.disc1Vec- mean(this.usVars.raw.disc1Vec); %[mm]

                    this.usVars.tVec = this.usVars.raw.tVec;
            end
        end

        function acCoupling(this)
           fprintf("US-Analysis: AC Coupling \n");
           this.us.dataAC =  this.us.data - mean(this.us.data, 4);
        end
        
        function extractP2P(this)
            fprintf("US-Analysis: Extractinf Spatial P2P \n");
            this.us.p2pMV = max(this.us.dataAC, [], 4) - min(this.us.dataAC, [], 4);
        end
        
        function convertToPascal(this)
            fprintf("US-Analysis: Converting To Pa \n");
            this.us.p2pRaw = this.us.p2pMV/837 *1e3;
        end
        
        function interpolate(this)
            fprintf("US-Analysis: Interpolating p2p \n");
            int = this.vars.intFactor;

            this.usVars.dx = this.usVars.dxRaw/int(1);
            this.usVars.dy = this.usVars.dyRaw/int(2);
            this.usVars.dz = this.usVars.dzRaw/int(3);
            
            this.usVars.xVec = this.usVars.xVecRaw(1) : this.usVars.dx : this.usVars.xVecRaw(end);
            this.usVars.yVec = this.usVars.yVecRaw(1) : this.usVars.dy : this.usVars.yVecRaw(end);
            this.usVars.zVec = this.usVars.zVecRaw(1) : this.usVars.dz : this.usVars.zVecRaw(end);
            
            this.us.p2p = interp3(this.usVars.zVecRaw, this.usVars.yVecRaw,this.usVars.xVecRaw,...
                                  this.us.p2pRaw,...
                                  this.usVars.zVec', this.usVars.yVec, this.usVars.xVec', 'cubic');
        end
        
        function extractFocusIdx(this)
           fprintf("US-Analysis: Finding Focus Index \n");
           % Find index of maximal P2P value.
           maxVal = max(this.us.p2pRaw(:));
           ind = find(this.us.p2pRaw == maxVal);
           [this.usVars.yMidIdxRaw, this.usVars.zMidIdxRaw, this.usVars.xMidIdxRaw] = ...
               ind2sub(size(this.us.p2pRaw), ind);
           
           this.usVars.maxPressure = max(this.us.p2p(:));
           ind = find(this.us.p2p == this.usVars.maxPressure);
           [this.usVars.yMidIdx, this.usVars.zMidIdx, this.usVars.xMidIdx] = ...
               ind2sub(size(this.us.p2p), ind);
        end

        function extractProfiles(this)
            fprintf("US-Analysis: Extracting Focal Profiles \n");
            this.us.trans1Profile =  squeeze(this.us.p2p(:, this.usVars.zMidIdx, this.usVars.xMidIdx));
            this.us.trans2Profile =  squeeze(this.us.p2p(this.usVars.yMidIdx, :, this.usVars.xMidIdx));
            this.us.depthProfile  =  squeeze(this.us.p2p(this.usVars.yMidIdx, this.usVars.zMidIdx, :));
            
            this.us.trans1ProfileNorm = analysisFunctions.normMatf(this.us.trans1Profile);
            this.us.trans2ProfileNorm = analysisFunctions.normMatf(this.us.trans2Profile);
            this.us.depthProfileNorm  = analysisFunctions.normMatf(this.us.depthProfile);
            
            this.us.trans1ProfileDB = db(this.us.trans1ProfileNorm);
            this.us.trans2ProfileDB = db(this.us.trans2ProfileNorm);
            this.us.depthProfileDB  = db(this.us.depthProfileNorm);
        end
        
        function extractSpeedOfSound(this)
            fprintf("US-Analysis: Estimating Speed-of-Sound \n");
            sigF = squeeze(this.us.data(this.usVars.yMidIdxRaw, this.usVars.zMidIdxRaw, :, :));
            for i=1:20
                [~, idxT(i)] = max(abs(sigF(this.usVars.xMidIdxRaw +(i-1),:)));
            end
            dTpeak =  abs(mean(this.usVars.tVec(idxT(2:end)) - this.usVars.tVec(idxT(1:end-1))));
            this.usVars.c = this.usVars.dxRaw*1e-3/dTpeak;
        end
        
        function calcDepthVec(this)
            fprintf("US-Analysis: Creating Depth Vec \n");
            this.usVars.depthVec = this.usVars.c * this.usVars.tVec *1e3; %[mm]
            this.usVars.dDepth = abs(this.usVars.depthVec(2) - this.usVars.depthVec(1)); 
        end
        
        function calcMinimalDelay(this)
            fprintf("US-Analysis: Finding Minimal Delay \n");
            sig = squeeze(this.us.dataAC(this.usVars.yMidIdxRaw, this.usVars.zMidIdxRaw, 1, :));
            sigG = gradient(sig);
            this.usVars.nearSigStartIdx = find(sigG(51:end) > 0.01*max(sigG),1) + 50;
            this.usVars.shortestDelay = this.usVars.tVec(this.usVars.nearSigStartIdx);
        end
        
        function fixUSAxisOffset(this)
            fprintf("US-Analysis: Fixing US axis (scan) \n");
            this.usVars.offset = this.usVars.c * this.usVars.shortestDelay *1e3; %[mm]
            this.usVars.xVecFix = this.usVars.xVec - min(this.usVars.xVec) + this.usVars.offset;
        end
        
        function calcFocalLen(this)
            fprintf("US-Analysis: Calculating Focal Length \n");
            this.usVars.focalLen = this.usVars.xVecFix(this.usVars.xMidIdx);
        end
        
        function calcRayleighLen(this)
            fprintf("US-Analysis: Calculating Rayleigh Length \n");
            this.usVars.rayleigh = sum(this.us.depthProfileNorm>=0.5) * this.usVars.dx;
        end
        
        function calcWaistSize(this)
            fprintf("US-Analysis: Calculating Waist Size \n");
            waist1Idx = sum(this.us.trans1ProfileNorm>=0.5);
            waist2Idx = sum(this.us.trans2ProfileNorm>=0.5);
            waist1 =  waist1Idx * this.usVars.dy;
            waist2 =  waist2Idx * this.usVars.dz;
            
            this.usVars.waist1FullIdx = waist1Idx;
            this.usVars.waist1Idx = ceil(waist1Idx/2);
            this.usVars.waist2FullIdx = waist2Idx;
            this.usVars.waist2Idx = ceil(waist2Idx/2);
            this.usVars.waistFull = (waist2 + waist1) / 2;
            this.usVars.waist = (waist2 + waist1) / 4;

        end
        
        function calcFocalFFT(this)
            fprintf("US-Analysis: FFT of Focus Temporal Profile \n");
            sig = squeeze(this.us.dataAC(this.usVars.yMidIdxRaw, this.usVars.zMidIdxRaw, this.usVars.xMidIdxRaw,:));
            this.us.focalFFT = abs(fftshift(fft(sig)));
            N = length(sig);
            this.usVars.fBar = (this.usVars.rawGenVars.fs/N) *  ( (-N/2) : 1 : (N/2)-1 );
        end
        
        function calcEnvelopes(this)
            fprintf("US-Analysis: Calculating Envelopes \n");
            if this.usVars.calcEnv
                envelopes = zeros(size(this.us.dataAC));
                curEnv = zeros(1, size(this.us.dataAC,2), size(this.us.dataAC,3), size(this.us.dataAC,4));
                for i=1:this.usVars.raw.scanSizeBin(1)
                    curData = this.us.dataAC(i,:,:,:);
                    parfor j=1:this.usVars.raw.scanSizeBin(2) 
                           curEnv(1,j,:,:) = permute(envelope(squeeze(curData(1,j,:,:))', 3, 'peak'), [3,4,2,1]); 
                    end
                    envelopes(i,:,:,:) = curEnv;
                end
                this.us.envelopes = envelopes;
                this.usVars.calcEnv = false;
            end
        end
        
        function extractFocalSignal(this)
            this.us.focalSig    = squeeze(this.us.dataAC(this.usVars.yMidIdxRaw, this.usVars.zMidIdxRaw, this.usVars.xMidIdxRaw, :));
            this.us.focalSigEnv = squeeze(this.us.envelopes(this.usVars.yMidIdxRaw, this.usVars.zMidIdxRaw, this.usVars.xMidIdxRaw,:));
        end
        
        
        
        function cutFocalData(this)
            fprintf("US-Analysis: Cutting Focal Data\n")
            this.us.focusEnvelope = this.us.envelopes((this.usVars.yMidIdx-this.usVars.waist1Idx) : (this.usVars.yMidIdx+this.usVars.waist1Idx),...
                                                      (this.usVars.zMidIdx-this.usVars.waist2Idx) : (this.usVars.zMidIdx+this.usVars.waist2Idx),...
                                                       :,:);
            this.us.focusPulses = this.us.dataAC((this.usVars.yMidIdx-this.usVars.waist1Idx) : (this.usVars.yMidIdx+this.usVars.waist1Idx),...
                                                      (this.usVars.zMidIdx-this.usVars.waist2Idx) : (this.usVars.zMidIdx+this.usVars.waist2Idx),...
                                                       :,:);
            this.usVars.yVecCut = (1:1:size(this.us.focusEnvelope,1))*this.usVars.dy;
            this.usVars.yVecCut = this.usVars.yVecCut - mean(this.usVars.yVecCut);
            this.usVars.zVecCut = (1:1:size(this.us.focusEnvelope,2))*this.usVars.dz;
            this.usVars.zVecCut = this.usVars.zVecCut - mean(this.usVars.zVecCut);
            
            this.usVars.yFocalIdxLen = length(this.usVars.yVecCut);
            this.usVars.zFocalIdxLen = length(this.usVars.zVecCut);
            this.usVars.yFocalIdx = ceil(this.usVars.yFocalIdxLen/2);
            this.usVars.zFocalIdx = ceil(this.usVars.zFocalIdxLen/2);
        end
        
        function interpolateFocalEnvelopes(this)
            fprintf("US-Analysis: Interpolating Focal Envelopes \n");
            this.us.focalEnvInt = [];
            this.us.focalPulseInt = [];
            
            focusEnvelope = this.us.focusEnvelope; 
            focusPulses   = this.us.focusPulses;
            
            % Match Axial resolution:
            fprintf("US-Analysis: Matching US axial resolution \n")
            this.usVars.dAx   = this.usVars.dDepth;
            this.usVars.axVec = this.usVars.xVecFix(1) : this.usVars.dAx : this.usVars.xVecFix(end);

            [Z, X, D] = meshgrid(this.usVars.xVecFix, this.usVars.zVecCut,  this.usVars.depthVec);
            [Zq, Xq, Dq] = meshgrid(this.usVars.axVec, this.usVars.zVecCut, this.usVars.depthVec);
            
            focalEnvIntAx = zeros(length(this.usVars.yVecCut),...
                           length(this.usVars.zVecCut),...
                           length(this.usVars.axVec),...
                           length(this.usVars.depthVec));
            focalPulseIntAX = zeros(length(this.usVars.yVecCut),...
                             length(this.usVars.zVecCut),...
                             length(this.usVars.axVec),...
                             length(this.usVars.depthVec));
            
            parfor i=1:length(this.usVars.yVecCut)
                focalEnvIntAx(i,:,:,:)   = interp3(Z, X, D, squeeze(focusEnvelope(i,:,:,:)), Zq, Xq, Dq, 'cubic');
                focalPulseIntAX(i,:,:,:) = interp3(Z, X, D, squeeze(focusPulses(i,:,:,:)), Zq, Xq, Dq, 'cubic');
            end

            [~, this.usVars.focalIndIntAx]    = min(abs(this.usVars.axVec - this.usVars.focalLen));
            [~, this.usVars.focalIndIntDepth] = min(abs(this.usVars.depthVec - this.usVars.focalLen));
            
            % Match transversal resolution:
            fprintf("US-Analysis: Matching US transversal resolution\n")
            if this.usVars.dy ~= this.usVars.dz
                if this.usVars.dy < this.usVars.dz
                    this.usVars.dtr1   = this.usVars.dy;
                    this.usVars.dtr2   = this.usVars.dy;
                    this.usVars.tr1Vec = this.usVars.yVecCut;
                    this.usVars.tr2Vec = this.usVars.zVecCut(1) : this.usVars.dtr2 : this.usVars.zVecCut(end);
                    
                    this.usVars.tr1FocalIdx = this.usVars.yFocalIdx;
                    [~, this.usVars.tr2FocalIdx] = min(abs( this.usVars.tr2Vec - this.usVars.zVecCut(this.usVars.zFocalIdx)));
                    
                    [Z, X, D] = meshgrid(this.usVars.axVec, this.usVars.zVecCut , this.usVars.depthVec);
                    [Zq, Xq, Dq] = meshgrid(this.usVars.axVec, this.usVars.tr2Vec, this.usVars.depthVec);
                    
                    focalEnvIntTR = zeros(length(this.usVars.tr1Vec),...
                                        length(this.usVars.tr2Vec),...
                                        length(this.usVars.axVec),...
                                        length(this.usVars.depthVec));
                    
                    focalPulseIntTR = zeros(length(this.usVars.tr1Vec),...
                                          length(this.usVars.tr2Vec),...
                                          length(this.usVars.axVec),...
                                          length(this.usVars.depthVec));
                    
                    parfor i=1:length(this.usVars.tr1Vec)
                        focalEnvIntTR(i,:,:,:) = interp3(Z, X, D, squeeze(focalEnvIntAx(i,:,:,:)), Zq, Xq, Dq, 'spline');
                        focalPulseIntTR(i,:,:,:) = interp3(Z, X, D,  squeeze(focalPulseIntAX(i,:,:,:)), Zq, Xq, Dq, 'spline');
                    end 
                    
%                      figure(); 
%                      ax(1) = subplot(1,2,1);
%                      imagesc(this.usVars.depthVec, this.usVars.xVecFix, squeeze(focusPulses(this.usVars.yFocalIdx,this.usVars.zFocalIdx,:,:)))
%                      title("PreInt")
%                      ax(2) = subplot(1,2,2);
%                      imagesc(this.usVars.depthVec, this.usVars.axVec, squeeze(focalPulseIntTR(this.usVars.tr1FocalIdx,this.usVars.tr2FocalIdx,:,:)))
%                      title("PostInt")
%                      linkaxes(ax)
                else
                    this.usVars.dtr2 = this.usVars.dz;
                    this.usVars.dtr1 = this.usVars.dz;
                    this.usVars.tr2Vec = this.usVars.zVecCut;
                    this.usVars.tr1Vec = this.usVars.yVecCut(1) : this.usVars.us.dtr1 : this.usVars.yVecCut(end);
                    
                    [~, this.usVars.tr1FocalIdx] = min(abs( this.usVars.tr1Vec - this.usVars.yVecCut(this.usVars.yFocalIdx)));
                    this.usVars.tr2FocalIdx = this.usVars.zFocalIdx;
                    
                    [Y, X, D] = meshgrid(this.usVars.axVec, this.usVars.yVecCut , this.usVars.depthVec);
                    [Yq, Xq, Dq] = meshgrid(this.usVars.axVec, this.usVars.tr1Vec, this.usVars.depthVec);
                    
                    focalEnvIntTR = zeros(length(this.usVars.tr1Vec),...
                                        length(this.usVars.tr2Vec),...
                                        length(this.usVars.axVec),...
                                        length(this.usVars.depthVec));
                    
                    focalPulseIntTR = zeros(length(this.usVars.tr1Vec),...
                                          length(this.usVars.tr2Vec),...
                                          length(this.usVars.axVec),...
                                          length(this.usVars.depthVec));

                    for i=1:length(this.usVars.tr2Vec)
                        focalEnvIntTR(:,i,:,:) = interp3(Y, X, D, squeeze(focalEnvIntAx(:,i,:,:)), Yq, Xq, Dq, 'cubic');
                        focalPulseIntTR(:,i,:,:) = interp3(Y, X, D, squeeze(focalPulseIntAX(:,i,:,:)), Yq, Xq, Dq, 'cubic');
                    end
                end
            end
                        
            this.us.focalEnvInt   = focalEnvIntTR;
            this.us.focalPulseInt = focalPulseIntTR;
        end
        
        function createSpatialData(this)
            fprintf("US-Analysis: Creating Spatial Data\n")
            this.us.spatialEnv           = permute(this.us.focalEnvInt, [1,2,4,3]);
            this.us.spatialPulse         = permute(this.us.focalPulseInt, [1,2,4,3]);
            this.usVars.spatial1Vec      = this.usVars.depthVec;
            this.usVars.spatial2Vec      = this.usVars.axVec;
            this.usVars.distFromTransVec = this.usVars.depthVec;
            
            [~, this.usVars.focalSpat1Idx] = min(abs(this.usVars.spatial1Vec - this.usVars.focalLen));
            [~, this.usVars.focalSpat2Idx] = min(abs(this.usVars.spatial2Vec - this.usVars.focalLen));
            
            figure()
            imagesc(this.usVars.spatial2Vec, this.usVars.spatial1Vec, squeeze(this.us.spatialEnv(9,9,:,:)))
        end
        
        function calcPulseWidth(this)
            fprintf("US-Analysis: Calculating Spatial Pulse Width\n")
            manualOffset = 0;
            for i=1:50
                curData = squeeze(this.us.spatialEnv(this.usVars.tr1FocalIdx,... 
                                                     this.usVars.tr2FocalIdx,...
                                                     :,...
                                                     this.usVars.focalSpat2Idx+i-1));
                maxVal = max(curData);
                pulseSize(i) = sum(curData > 0.01*maxVal);
            end
            this.usVars.pulseSizeIdx = ceil(mean(pulseSize))+manualOffset;
            this.usVars.pulseSizeLen = this.usVars.pulseSizeIdx*this.usVars.dAx;
        end

        function cutPulses(this)
            fprintf("US-Analysis: Cutting Pulses \n")
            % 1: crop only data in relevant region:
            startIdx = find(this.usVars.spatial1Vec == this.usVars.spatial2Vec(1));
            endIdx   = find(this.usVars.spatial1Vec == this.usVars.spatial2Vec(end));
            
            this.us.focalEnvCrop     = this.us.spatialEnv(:,:,startIdx:endIdx, :); 
            this.us.focalPulseCrop   = this.us.spatialPulse(:,:,startIdx:endIdx,:); 
            
            this.usVars.spatial1VecCrop = this.usVars.spatial1Vec(startIdx:endIdx);
            % this.usVars.spatial1VecCrop = this.usVars.spatial1Vec(startIdx+this.usVars.pulseSizeIdx:endIdx);
            
            this.usVars.numOfCutPulses = length(this.usVars.spatial1VecCrop)-this.usVars.pulseSizeIdx+1;
            
            % 2: Cut-out pulses from cropped data:
            this.us.focalPulsesEnvCut = zeros(length(this.usVars.tr1Vec),...
                                              length(this.usVars.tr2Vec),...
                                              this.usVars.numOfCutPulses,...
                                              this.usVars.pulseSizeIdx);
            this.us.focalPulsesCut    = zeros(length(this.usVars.tr1Vec),...
                                              length(this.usVars.tr2Vec),...
                                              this.usVars.numOfCutPulses,...
                                              this.usVars.pulseSizeIdx);
            
            this.usVars.pulseVec = this.usVars.spatial1Vec(1:this.usVars.pulseSizeIdx);
            this.usVars.pulseStartPosVec = this.usVars.spatial1VecCrop(this.usVars.pulseSizeIdx:end);
                                          
            for i=1:this.usVars.numOfCutPulses
                idxStart= i + this.usVars.pulseSizeIdx -1;
                idxStartLine = idxStart - this.usVars.pulseSizeIdx +1;
                this.us.focalPulsesEnvCut(:,:,i,:) = this.us.focalEnvCrop(:,:,idxStart,idxStartLine:idxStart);
                this.us.focalPulsesCut(:,:,i,:)    = this.us.focalPulseCrop(:,:,idxStart,idxStartLine:idxStart);
            end
%             figure();
%             subplot(1,2,1)
%             plot(squeeze( this.us.focalEnvCrop(9,9,i+this.usVars.pulseSizeIdx-1,:)));
%             subplot(1,2,2)
%             plot(squeeze(this.us.focalPulsesEnvCut(9,9,i,:)))
        end
        
        function calaPulsesProfile(this)
            
            this.us.focalProfile = max(this.us.focalPulsesEnvCut,[],4);

%             figure()
%             imagesc(squeeze(this.us.focalProfile(9,:,:)))
            
        end
        
        function displayResults(this)
            fprintf("US-Analysis: Displaying Results\n");
            %------------------------------
            %1D profiles:
            %------------------------------
            figure();
            subplot(3,2,1)
            plot(this.usVars.yVec, this.us.trans1ProfileNorm); hold on
            plot(this.usVars.yVec, 0.5*ones(length(this.usVars.yVec),1))
            xlabel("TR1[mm]")
            ylabel("Lin. Power [AU] ")
            subplot(3,2,3)
            plot(this.usVars.zVec, this.us.trans2ProfileNorm);hold on
            plot(this.usVars.zVec, 0.5*ones(length(this.usVars.zVec),1))
            xlabel("TR2[mm]")
            ylabel("Lin. Power [AU] ")
            subplot(3,2,5)
            plot(this.usVars.xVecFix, this.us.depthProfileNorm);hold on
            plot(this.usVars.xVecFix, 0.5*ones(length(this.usVars.xVecFix),1))
            xlabel("Depth [mm]")
            ylabel("Lin. Power [AU] ")
            subplot(3,2,2)
            plot(this.usVars.yVec, this.us.trans1ProfileDB);hold on
            plot(this.usVars.yVec, -6*ones(length(this.usVars.yVec),1))
            xlabel("TR1[mm]")
            ylabel("Log. Power [dB]")
            subplot(3,2,4)
            plot(this.usVars.zVec, this.us.trans2ProfileDB);hold on
            plot(this.usVars.zVec, -6*ones(length(this.usVars.zVec),1))
            xlabel("TR2[mm]")
            ylabel("Log. Power [dB]")
            subplot(3,2,6)
            plot(this.usVars.xVecFix, this.us.depthProfileDB);hold on
            plot(this.usVars.xVecFix, -6*ones(length(this.usVars.xVecFix),1))
            xlabel("Depth [mm]")
            ylabel("Log. Power [dB]")
            
            %------------------------------
            %2D Profiles:
            %------------------------------
            figure(); 
            subplot(3,1,1)
            imagesc(this.usVars.xVecFix,...
                    this.usVars.yVec,...
                    squeeze(this.us.p2p(:,this.usVars.zMidIdx,:)))
            axis tight equal
            xlabel("US[mm]")
            ylabel("Trans 1 [mm]")
            hCB = colorbar;
            ylabel(hCB, "KPa")
            subplot(3,1,2)
            imagesc(this.usVars.xVecFix,...
                    this.usVars.zVec,...
                    squeeze(this.us.p2p(this.usVars.yMidIdx,:,:)))
            xlabel("US[mm]")
            ylabel("Trans 2 [mm]")
            axis tight equal
            hCB = colorbar;
            ylabel(hCB, "KPa")
            subplot(3,1,3)
            imagesc(this.usVars.yVec,...
                    this.usVars.zVec,...
                    squeeze(this.us.p2p(:,:,this.usVars.xMidIdx)))
            xlabel("Trans 1[mm]")
            ylabel("Trans 2 [mm]")
            axis tight equal
            hCB = colorbar;
            ylabel(hCB, "KPa")
            
            %------------------------------
            % Focal Point Pulse: Temporal and FFT
            %------------------------------
            figure()
            subplot(1,2,1)
            plot(this.usVars.tVec*1e6, squeeze(this.us.dataAC(this.usVars.yMidIdxRaw, this.usVars.zMidIdxRaw, this.usVars.xMidIdxRaw,:))); hold on
            plot(this.usVars.tVec*1e6, squeeze(this.us.envelopes(this.usVars.yMidIdxRaw, this.usVars.zMidIdxRaw, this.usVars.xMidIdxRaw,:)));
            xlabel("t[\mu s]");
            ylabel("Amplitude [AU]");
            subplot(1,2,2)
            plot(this.usVars.fBar*1e-6, this.us.focalFFT);
            xlabel("f[MHz]");
            ylabel("Amplitude [AU]");
            
            %------------------------------
            % Focal Point Pulse: Temporal and FFT
            %------------------------------
            figure()
            NIm = this.usVars.raw.scanSizeBin(3);
            dIm = floor(NIm/10);
            for i=1:dIm:NIm
                plot(this.usVars.tVec*1e6, squeeze(this.us.dataAC(this.usVars.yMidIdxRaw, this.usVars.zMidIdxRaw, i,:))); hold on
                plot(this.usVars.tVec*1e6, squeeze(this.us.envelopes(this.usVars.yMidIdxRaw, this.usVars.zMidIdxRaw, i,:))); hold on
            end
            xlabel("t[\mu s]");
            ylabel("Amplitude [AU]");
            title("Pulses And Envelopes along Focal Axis (Pre-Interpolation)")

            
            %------------------------------
            % 2D Pulses and Envelopes along Focal Axis (Pre Interpolation):
            %------------------------------
            figure();
            subplot(1,2,1)
            imagesc(this.usVars.tVec*1e6, this.usVars.xVecFix,...
                    squeeze(this.us.dataAC(this.usVars.yMidIdxRaw, this.usVars.zMidIdxRaw, :,:)));
            colorbar;
            axis tight equal
            xlabel("t[\mu s]");
            ylabel("X [mm]");
            title("Pulses")
            subplot(1,2,2)
            imagesc(this.usVars.tVec*1e6, this.usVars.xVecFix,...
                    squeeze(this.us.envelopes(this.usVars.yMidIdxRaw, this.usVars.zMidIdxRaw, :,:)));
            colorbar;
            axis tight equal
            xlabel("t[\mu s]");
            ylabel("X [mm]");
            title("Envelopes")
            
            %------------------------------
            % Print Beam Parameters:
            %------------------------------
            Parameter = ["Focal Length"; "Rayleigh Length"; "Waist Size"; "Maximun Pressure"];
            Value     = [this.usVars.focalLen; this.usVars.rayleigh; this.usVars.waist; this.usVars.maxPressure];
            Units = ["mm" ; "mm"; "mm"; "KPa"];
            T = table(Parameter, Value, Units);
            disp(T);
            
            %------------------------------
            % Interpolation
            %------------------------------
            figure()
            
            ax(1) = subplot(2,2,1);
            imagesc(this.usVars.depthVec, this.usVars.xVecFix,...
                    squeeze(this.us.dataAC(this.usVars.yMidIdxRaw, this.usVars.zMidIdxRaw, :,:)));
            colorbar;
            axis tight equal
            xlabel("Depth[mm]");
            ylabel("X [mm]");
            title("Pulses - Pre Int.")
            
            ax(2) = subplot(2,2,2);
            imagesc(this.usVars.depthVec, this.usVars.xVecFix,...
                    squeeze(this.us.envelopes(this.usVars.yMidIdxRaw, this.usVars.zMidIdxRaw, :,:)));
            colorbar;
            axis tight equal
            xlabel("Depth[mm]");
            ylabel("X [mm]");
            title("Envelopes - Pre Int.")
            
            size1 = this.usVars.tr1FocalIdx;
            size2 = this.usVars.tr2FocalIdx;
            
            ax(3) = subplot(2,2,3);
            imagesc(this.usVars.depthVec, this.usVars.axVec, ...
                    squeeze(this.us.focalPulseInt(size1, size2, :,:)));
            colorbar;
            axis tight equal
            xlabel("Depth[mm]");
            ylabel("X [mm]");
            title("Pulses - Post Int.")
            
            ax(4) = subplot(2,2,4);
            imagesc(this.usVars.depthVec, this.usVars.axVec, ...
                    squeeze(this.us.focalEnvInt(size1, size2, :,:)));
            colorbar;
            axis tight equal
            xlabel("Depth[mm]");
            ylabel("X [mm]");
            title("Envelopes - Post Int.")
            
            linkaxes(ax);
            
            %------------------------------
            %Spatial Data (1D)          
            %------------------------------
            figure()
            size1 = this.usVars.tr1FocalIdx;
            size2 = this.usVars.tr2FocalIdx;
            subplot(1,2,1)
            for i =1:100:length(this.usVars.axVec)
                plot(this.usVars.tVec*1e6, squeeze(this.us.focalEnvInt(size1,size2,i,:))); hold on
                plot(this.usVars.tVec*1e6, squeeze(this.us.focalPulseInt(size1,size2,i,:))); 
            end
            xlabel ("Time [\mu s]")
            hold off;

            subplot(1,2,2)
            for i =1:100:length(this.usVars.spatial1Vec)
%                 if i==1; i=20; end
                plot(this.usVars.spatial2Vec, squeeze(this.us.spatialEnv(size1,size2,i,:))); hold on
                plot(this.usVars.spatial2Vec, squeeze(this.us.spatialPulse(size1,size2,i,:)));
            end
            hold off; 
            xlabel ("X [mm]")
            
            
            %------------------------------
            %Interpolation & Spatial Data (2D)          
            %------------------------------
            figure()
            subplot(2,2,1)
            imagesc(this.usVars.depthVec, this.usVars.axVec, squeeze(this.us.focalEnvInt(size1,size2,:,:)))
            axis tight equal
            title( "Envelope Sinogram")
            xlabel ("Distance From Trans. [mm] [mm]")
            ylabel ("Depth [mm]")
            
            subplot(2,2,2)
            imagesc(this.usVars.spatial2Vec, this.usVars.spatial2Vec, squeeze(this.us.spatialEnv(size1,size2,:,:)))
            axis tight equal
            title( "Envelope Spatial")
            xlabel ("Depth [mm]")
            ylabel ("Distance From Trans. [mm]")
            colorbar
            
            subplot(2,2,3)
            imagesc(this.usVars.depthVec, this.usVars.axVec, squeeze(this.us.focalPulseInt(size1,size2,:,:)))
            axis tight equal
            title( "Pulses Sinogram")
            xlabel ("Distance From Trans. [mm] [mm]")
            ylabel ("Depth [mm]")
            
            subplot(2,2,4)
            imagesc(this.usVars.spatial2Vec, this.usVars.spatial2Vec, squeeze(this.us.spatialPulse(size1,size2,:,:)))
            axis tight equal
            title( "Pulses Spatial")
            xlabel ("Depth [mm]")
            ylabel ("Distance From Trans. [mm]")
            colorbar
            
            %------------------------------
            %Interpolation & Spatial Data           
            %------------------------------
%             lenY = size(this.us.focusEnvelope,1);
%             lenZ = size(this.us.focusEnvelope,2);
%             figure()
%             subplot(2,2,1)
%             for i =1:100:length(this.usVars.depthAxPulses)
%                 if i==1; i=20; end
%                 plot(squeeze(this.us.pulses(size1,size2,i,:))); hold on
%             end
%             
%             subplot(2,2,2)
%             for i =1:2:lenY
%                 plot(squeeze(this.us.pulses(i,size2,370,:))); hold on
%             end
%             
%             subplot(2,2,3)
%             for i =1:2:lenZ
%                 plot(squeeze(this.us.pulses(size1,i,370,:))); hold on
%             end
%             xlabel(
        end
        
    end
end