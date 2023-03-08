classdef usLoader < handle
    %USLOADER Summary of this class goes hereLoader
    %   Detailed explanation goes here
    % Separate measurements of depth profile and 2D focal transverse
    % profile.
    % **Assumption: the temporal sampling discretization is identical.**
    
    properties
        uVars
        data
        res
        scanVars
        vars
        grid
        usStats
    end
    
    methods (Static)
        function uVars = createUserVars()
            uVars.paths = {};
            uVars.intFactor      = [];
            uVars.usAtZero       = false;
            uVars.useExtSpeed    = false;
            uVars.extC           = [];
            uVars.pulseSizeIdx   = 100;
            uVars.envPeakSize    = 10;
            uVars.displayResults = false;
        end
    end

    methods
        function this = usLoader()

        end

        function res = loadAndAnalyse(this, uVars)
            this.uVars = uVars;
            this.loadDB();

            this.analyse();

            res.data = this.data;
            res.grid = this.grid;
            res.usStats = this.usStats;

            if this.uVars.displayResults
                this.displayResults();
            end
        end

        function loadDB(this)
            %load Depth scan
            tmp = load(this.uVars.paths{1});
            close all;
            this.data.raw.depth = squeeze(tmp.resCs);
            if ~this.uVars.usAtZero
                this.data.raw.depth = flip(this.data.raw.depth,1);
            end

            this.scanVars.depth = tmp.csVars;
            
            this.grid.raw.axialVecRaw = this.scanVars.depth.scanVecBin;
            this.grid.raw.dAx      = this.scanVars.depth.axScanStride;
            this.grid.raw.tVec     = this.scanVars.depth.tVec;
            this.grid.raw.dt       = 1/this.scanVars.depth.fs;
            this.vars.fs           = this.scanVars.depth.fs;

            %load transvers scan
            tmp = load(this.uVars.paths{2});
            close all;
            this.data.raw.tr = squeeze(tmp.resCs);
            this.scanVars.tr = tmp.csVars;

            this.grid.raw.tr1VecRaw = this.scanVars.tr.scanVecBin;
            this.grid.raw.tr2VecRaw = this.scanVars.tr.disc1Vec;
            this.grid.raw.dtr1   = this.scanVars.tr.axScanStride;
            this.grid.raw.dtr2   = this.scanVars.tr.axDisc1Stride;
            
            this.grid.raw.axialSize = length(this.grid.raw.axialVecRaw);
            this.grid.raw.tr1Size   = length(this.grid.raw.tr1VecRaw);
            this.grid.raw.tr2Size   = length(this.grid.raw.tr2VecRaw);

            % Collect commin variables:
            this.grid.tVec             = this.grid.raw.tVec;
            this.grid.dt               = this.grid.raw.dt;
        end
        
        function analyse(this)
            this.acCoupling();
            this.convertToPascal();
            this.extractP2P();
            this.extractFocusIdx();
            this.extractSpeedOfSound();

            this.interpolate();
            this.alignAxialAndDepthVec();
            this.calcEnvelopes();

            this.calcBeamGeometry();
            this.extractFocalPulse();

            this.createSpatialData();
            this.cutPulses();
        end

        function acCoupling(this)
            fprintf("US-Loader: AC Coupling & mV conversion \n");
            this.data.raw.depthAC = (this.data.raw.depth - mean(this.data.raw.depth, 2) ) * 1e3;
            this.data.raw.trAC    = (this.data.raw.tr    - mean(this.data.raw.tr, 3)) * 1e3;
        end
        
        function convertToPascal(this)
            fprintf("US-Loader: Converting To KPa \n");
            this.data.raw.depthPa = this.data.raw.depthAC/837 *1e3; %KPa
            this.data.raw.trPa    = this.data.raw.trAC/837 *1e3;    %KPa
        end

        function extractP2P(this)
            fprintf("US-Loader: Extracting Spatial P2P \n");
            this.data.raw.depthP2p = max(this.data.raw.depthPa, [], 2) - min(this.data.raw.depthPa, [], 2);
            this.data.raw.trP2p    = max(this.data.raw.trPa, [], 3)    - min(this.data.raw.trPa, [], 3);
        end

        function extractFocusIdx(this)
            fprintf("US-Loader: Finding Focus Index \n");
            % Find index of maximal P2P value:
            [maxValAxial, axialFocalIdx] = max(this.data.raw.depthP2p);
            
            maxValTr = max(this.data.raw.trP2p(:));
            
            ind = find(this.data.raw.trP2p == maxValTr);
            [tr1Ind, tr2Ind] = ind2sub(size(this.data.raw.trP2p), ind);
            
            this.usStats.maxPressureRaw = max([maxValAxial, maxValTr]);
            
            this.grid.raw.axialFocalIdx = axialFocalIdx;
            this.grid.raw.tr1FocalIdx   = tr1Ind;
            this.grid.raw.tr2FocalIdx   = tr2Ind;

            this.grid.raw.tr1Vec = this.grid.raw.tr1VecRaw - this.grid.raw.tr1VecRaw(tr1Ind);
            this.grid.raw.tr2Vec = this.grid.raw.tr2VecRaw - this.grid.raw.tr2VecRaw(tr2Ind);
        end

        function extractSpeedOfSound(this)
            fprintf("US-Loader: Estimating Speed-of-Sound \n");
            depthSig   = this.data.raw.depthAC;
            dDpeth     = this.grid.raw.dAx;
            axFocalIdx = this.grid.raw.axialFocalIdx;
            tVec       = this.grid.tVec;

            calcLen = 100;
            idxT = zeros(1,calcLen);
            for i=1:calcLen
                [~, idxT(i)] = max(abs(depthSig(axFocalIdx +(i-1),:)));
            end
            dtPeakVec = tVec(idxT(2:end)) - tVec(idxT(1:end-1));
            dtPeak    = mean(dtPeakVec);
            this.usStats.cCalc = dDpeth*1e-3/dtPeak;

            if this.uVars.useExtSpeed
                this.usStats.c = this.uVars.extC;
            else
                this.usStats.c = this.usStats.cCalc;
            end

            this.grid.raw.depthVec  = this.usStats.c * this.grid.tVec * 1e3; %[mm]
            this.grid.raw.dDepth    = this.usStats.c * this.grid.dt   * 1e3; %[mm]
            this.grid.raw.depthSize = length(this.grid.raw.depthVec);
        end
       
        function interpolate(this)
            fprintf("US-Loader: Interpolating \n");
            int = this.uVars.intFactor;

            dAxInt  = this.grid.raw.dDepth; % match temporal and spatial resolution on scan and depth axes.
            dtr1Int = this.grid.raw.dtr1 / int(1);
            dtr2Int = this.grid.raw.dtr2 / int(2);
    
            this.grid.int.dAx  = dAxInt;
            this.grid.int.dtr1 = dtr1Int;
            this.grid.int.dtr2 = dtr2Int;
            this.grid.int.dDepth = this.grid.raw.dDepth;

            this.grid.int.axialVecRaw = this.grid.raw.axialVecRaw(1)  : dAxInt  : this.grid.raw.axialVecRaw(end);
            this.grid.int.tr1Vec      = this.grid.raw.tr1Vec(1)    : dtr1Int : this.grid.raw.tr1Vec(end);
            this.grid.int.tr2Vec      = this.grid.raw.tr2Vec(1)    : dtr2Int : this.grid.raw.tr2Vec(end);
            this.grid.int.depthVec    = this.grid.raw.depthVec;
            
            this.grid.int.axialSize = length(this.grid.int.axialVecRaw);
            this.grid.int.tr1Size   = length(this.grid.int.tr1Vec);
            this.grid.int.tr2Size   = length(this.grid.int.tr2Vec);
            this.grid.int.depthSize = length(this.grid.int.depthVec);

            % Interpolate Depth (sptail - temporal):
            [D, A]       = meshgrid(this.grid.raw.depthVec, this.grid.raw.axialVecRaw);
            [Dint, Aint] = meshgrid(this.grid.int.depthVec, this.grid.int.axialVecRaw);

            this.data.int.depth  = interp2( D, A, ...
                                            this.data.raw.depthPa,...
                                            Dint, Aint, 'cubic');

            this.data.int.depthP2p    = max(this.data.int.depth, [], 2) - min(this.data.int.depth, [], 2);
            this.data.int.depthP2pSmt = envelope(this.data.int.depthP2p, 10, 'peak'); 
            % Interpolate Transverse (Only P2P) :
            [TR2, TR1]       = meshgrid(this.grid.raw.tr1Vec, this.grid.raw.tr2Vec);
            [TR2Int, TR1Int] = meshgrid(this.grid.int.tr1Vec, this.grid.int.tr2Vec);

            this.data.int.trP2p  = interp2( TR2, TR1, ...
                                              this.data.raw.trP2p,...
                                              TR2Int, TR1Int,'cubic');

            this.data.int.trP2pNorm = analysisFunctions.normMatf(this.data.int.trP2p);
            this.data.int.trP2pDb   = db(this.data.int.trP2pNorm);

            % Find index of maximal *interpolated* P2P value:
            [maxValAxial, axialFocalIdx] = max(this.data.int.depthP2p);
            
            maxValTr = max(this.data.int.trP2p(:));
            
            ind = find(this.data.int.trP2p == maxValTr);
            [tr1Ind, tr2Ind] = ind2sub(size(this.data.int.trP2p), ind);
            
            this.usStats.maxPressure = max([maxValAxial, maxValTr]);
            
            this.grid.int.axialFocalIdx = axialFocalIdx;
            this.grid.int.tr1FocalIdx   = tr1Ind;
            this.grid.int.tr2FocalIdx   = tr2Ind;
            
            % Extract Interpolated profiles on main Axes:
            this.data.int.depthProf     = this.data.int.depthP2pSmt;
            this.data.int.depthProfNorm = this.data.int.depthProf/max(this.data.int.depthProf);
            this.data.int.depthProfDb   = db(normMatf(this.data.int.depthProf));

            this.data.int.tr1Prof     = this.data.int.trP2p(:, this.grid.int.tr2FocalIdx);
            this.data.int.tr1ProfNorm = this.data.int.trP2pNorm(:, this.grid.int.tr2FocalIdx);
            this.data.int.tr1ProfDb   = this.data.int.trP2pDb(:, this.grid.int.tr2FocalIdx);

            this.data.int.tr2Prof     = this.data.int.trP2p(this.grid.int.tr1FocalIdx, :);
            this.data.int.tr2ProfNorm = this.data.int.trP2pNorm(this.grid.int.tr1FocalIdx, :);
            this.data.int.tr2ProfDb   = this.data.int.trP2pDb(this.grid.int.tr1FocalIdx, :);

%             figure();
%             hp = plot(this.grid.int.depthVec, zeros(1,size(this.data.int.depth,2)));
%             ylim( [min(this.data.int.depth(:)), max(this.data.int.depth(:))] )
%             for i=1:size(this.data.int.depth, 1)
%                 set(hp, 'YData', this.data.int.depth(i, :));
%                 drawnow();
%             end
        end
        
        function alignAxialAndDepthVec(this)
            fprintf("US-Loader: Aligning Scan Depth Grid \n");

            axialVecRaw   = this.grid.raw.axialVecRaw;
            tVec          = this.grid.tVec;
            axialFocalIdx = this.grid.raw.axialFocalIdx;
            
            sig    = this.data.raw.depthPa(axialFocalIdx,:);
            [~, I] = max(abs(sig));
            
            envLow   = 40; envHigh = 100;
            envIdxs  = I-envLow:I+envHigh;
            pulseEnv = sig(envIdxs);
            envG     = gradient(pulseEnv);

            delayIdx = envIdxs(find(envG > 1 ,1));
            delay    = tVec(delayIdx);

            trueAxDepth = delay * this.usStats.c * 1e3;
            
            axialOffset = axialVecRaw(axialFocalIdx) - trueAxDepth;

            this.grid.raw.axialVec       = axialVecRaw - axialOffset;
            this.grid.raw.focalDelayIdx  = delayIdx;

            this.usStats.focalDelayRaw   = delay;
            

            axialVecRaw   = this.grid.int.axialVecRaw;
            tVec          = this.grid.tVec;
            axialFocalIdx = this.grid.int.axialFocalIdx;
            
            sig    = this.data.int.depth(axialFocalIdx,:);
            [~, I] = max(abs(sig));
            
            envLow   = 40; envHigh = 100;
            envIdxs  = I-envLow:I+envHigh;
            pulseEnv = sig(envIdxs);
            envG     = gradient(pulseEnv);

            delayIdx = envIdxs(find(envG > 1 ,1));
            delay    = tVec(delayIdx);

            trueAxDepth = delay * this.usStats.c * 1e3;
            
            axialOffset = axialVecRaw(axialFocalIdx) - trueAxDepth;

            this.grid.int.axialVec       = axialVecRaw - axialOffset;
            this.grid.int.focalDelayIdx  = delayIdx;

            this.usStats.focalDelay      = delay;
%             figure();
%             hp = plot(this.grid.raw.depthVec, zeros(1,size(this.data.raw.depthPa,2)));
%             ylim( [min(this.data.raw.depthPa(:)), max(this.data.raw.depthPa(:))] )
%             for i=1:1500
%                 set(hp, 'YData', this.data.raw.depthPa(i, :));
%                 drawnow();
%             end

        end

        function calcBeamGeometry(this)
            fprintf("US-Loader: Calculating Focal Length \n");
            this.usStats.focalLen = this.grid.int.axialVec(this.grid.int.axialFocalIdx);
            
            fprintf("US-Loader: Calculating Rayleigh Length \n");
            this.usStats.rayleigh = sum(this.data.int.depthProfNorm>=0.5) * this.grid.int.dAx;
            
            fprintf("US-Loader: Calculating Waist Size \n");
            waist1Idx = sum(this.data.int.tr1ProfNorm>=0.5);
            waist2Idx = sum(this.data.int.tr2ProfNorm>=0.5);
            waist1 =  waist1Idx * this.grid.int.dtr1;
            waist2 =  waist2Idx * this.grid.int.dtr2;

            this.usStats.waistFull = (waist2 + waist1) / 2;
            this.usStats.waist     = (waist2 + waist1) / 4;
            
            varsNames = {'Parameter' ; 'Values'; 'Units'};
            units =  ["mm" ; "mm"; "mm"; "KPa"; "m/s"];
            parNames  = ["Focal Length" ; "Waist Size"; "Rayleigh"; "Max Pressure"; "C"];
            values = [this.usStats.focalLen; this.usStats.waist; this.usStats.rayleigh; this.usStats.maxPressureRaw; this.usStats.c];
            this.usStats.T = table(parNames, values, units, 'VariableNames', varsNames); 
        end
        
        function calcEnvelopes(this)
            fprintf("US-Loader: Calculating Envelopes \n");
            this.data.int.depthEnv = envelope(this.data.int.depth', this.uVars.envPeakSize, 'peak'); 
            this.data.int.depthEnv =  this.data.int.depthEnv';
        end

        function extractFocalPulse(this)
            fprintf("US-Loader: FFT of Focus Temporal Profile \n");
            this.data.focal.pulse    = this.data.int.depth(this.grid.int.axialFocalIdx,:);
            this.data.focal.env      = this.data.int.depthEnv(this.grid.int.axialFocalIdx,:);
            this.grid.focal.depthVec = this.grid.int.depthVec;
            this.grid.focal.dDepth   = this.grid.int.dDepth;

            this.data.focal.FFT = abs(fftshift(fft(this.data.focal.pulse))).^2;
            this.data.focal.FFTNorm = this.data.focal.FFT/max(this.data.focal.FFT);
            N = this.grid.int.depthSize;
            fs = this.vars.fs;
            this.grid.focal.fBar = (fs/N) *  ( (-N/2) : 1 : (N/2)-1 );
        end
        
        function createSpatialData(this)
            depth    = this.data.int.depth;
            depthEnv = this.data.int.depthEnv;
            
            depthVec = this.grid.int.depthVec;
            axialVec = this.grid.int.axialVec;
            
            focalDelayIdx = this.grid.int.focalDelayIdx;
            axialFocalIdx = this.grid.int.axialFocalIdx;
            
            depthOffset = focalDelayIdx - axialFocalIdx;
            pulseSize   = 120;
            this.vars.pulseSizeRaw = pulseSize;
            
            axialI = 1;
            axialF = this.grid.int.axialSize;
            
            depthI = depthOffset + pulseSize;
            depthF = depthI      + this.grid.int.axialSize-pulseSize;
            
            this.data.spat.depth     = depth(axialI:axialF, depthI:depthF);
            this.data.spat.depthEnv  = depthEnv(axialI:axialF, depthI:depthF);
            
            this.grid.spat.axialVec  = axialVec(axialI:axialF);
            this.grid.spat.depthVec  = depthVec(depthI:depthF);
            
            this.grid.spat.depthSize = length(this.grid.spat.depthVec);
            this.grid.spat.axialSize = length(this.grid.spat.depthVec);
            
            this.data.spat.depthP2p = peak2peak(this.data.spat.depth, 2);
            this.grid.spat.axialFocalIdx = this.grid.int.axialFocalIdx;
        end

        function cutPulses(this)
            pulseSize = this.vars.pulseSizeRaw;

            pulses    = zeros(this.grid.spat.depthSize, pulseSize);
            pulsesEnv = zeros(this.grid.spat.depthSize, pulseSize);
            pulsesAx  = zeros(this.grid.spat.depthSize, pulseSize);

            for i=1:this.grid.spat.depthSize
                pulses(i,:)    = this.data.spat.depth(i:i+pulseSize-1, i);
                pulsesEnv(i,:) = this.data.spat.depthEnv(i:i+pulseSize-1, i);
                pulsesAx (i,:) = this.grid.int.axialVec(i:i+pulseSize-1);
            end
            
            axialFocalIdx = this.grid.spat.axialFocalIdx;

            this.data.pulses.pulses            = pulses;
            this.data.pulses.pulsesEnv         = pulsesEnv;
            this.data.pulses.focalPulse        = pulses(axialFocalIdx, :);
            this.data.pulses.focalPulseNorm    = this.data.pulses.focalPulse ./ max(abs(this.data.pulses.focalPulse));
            this.data.pulses.focalPulseEnv     = pulsesEnv(axialFocalIdx, :);
            this.data.pulses.focalPulseEnvNorm = this.data.pulses.focalPulseEnv / max(abs(this.data.pulses.focalPulseEnv));

            this.grid.pulses.axialVec      = this.grid.int.depthVec(1:pulseSize);
            this.grid.pulses.pulseStartPos = this.grid.spat.axialVec(pulseSize:end);
            this.grid.pulses.pulsesAx      = pulsesAx;
            this.grid.pulses.depthVec      = this.grid.spat.depthVec;

            this.grid.pulses.axialFocalIdx = axialFocalIdx;
            
%             figure();
%             ax1 = subplot(1,2,1);
%             hp1 = plot(ax1, this.grid.pulses.axialVec, zeros(1,pulseSize));
%             ylim( [min(pulses(:)), max(pulses(:))] )
%             ax2 = subplot(1,2,2);
%             hp2 = plot(ax2, this.grid.pulses.axialVec, zeros(1,pulseSize));
%             ylim( [min(pulsesEnv(:)), max(pulsesEnv(:))] )
%             for i=1:size(this.data.pulses.pulses, 1)
%                 set(hp1, 'YData', pulses(i,:));
%                 set(hp2, 'YData', pulsesEnv(i,:));
%                 drawnow();
%             end

        end

        function displayResults(this)
            axialVec = this.grid.int.axialVec;
            tr1 = this.grid.int.tr1Vec;
            tr2 = this.grid.int.tr2Vec;
            depthVec = this.grid.int.depthVec;
            fBar = this.grid.focal.fBar;
            
            disp(this.usStats.T);

            figure(); 
            subplot(3,2,1)
            plot(tr1, this.data.int.tr1Prof);
            xlabel("Transversal Axis 1 [mm]");
            ylabel("Pressure [KPa]")
            subplot(3,2,2)
            plot(tr1, this.data.int.tr1ProfDb);
            xlabel("Transversal Axis 1 [mm]");
            ylabel("Acoustic Intensity[Db]")
            subplot(3,2,3)
            plot(tr2,this.data.int.tr2Prof);
            xlabel("Transversal Axis 2 [mm]");
            ylabel("Pressure [KPa]")
            subplot(3,2,4)
            plot(tr2, this.data.int.tr2ProfDb);
            xlabel("Transversal Axis 2 [mm]");
            ylabel("Acoustic Intensity[Db]");
            subplot(3,1,3)
            imagesc(tr2, tr1, this.data.int.trP2p);
            axis tight equal
            xlabel("Transversal Axis 2 [mm]");
            ylabel("Transversal Axis 1 [mm]");
            hCB = colorbar();
            ylabel(hCB, 'Pressure [KPa]');
            
            figure();
            subplot(1,2,1)
            plot(axialVec, this.data.int.depthProf);
            xlabel("Depth Scan Axis [mm]")
            ylabel("Pressure [KPa]");
            subplot(1,2,2)
            plot(axialVec, this.data.int.depthProfDb);
            xlabel("Depth Scan Axis [mm]")
            ylabel("Acoustic Intensity[Db]");
            
            figure();
            subplot(1,2,1)
            plot(depthVec, this.data.focal.pulse); hold on
            plot(depthVec, this.data.focal.env);
            xlabel("Depth Scan Axis [mm]")
            ylabel("Pressure [KPa]");
            subplot(1,2,2)
            plot(fBar*1e-6, this.data.focal.FFTNorm)
            xlabel("Frequency [MHz]")
            ylabel("Acoustic Power Spectrum [AU]");
            

            figure();
            ax = subplot(2,2,1);
            imagesc(ax, this.grid.int.axialVec, this.grid.int.depthVec, this.data.int.depth);
%             imagesc(ax, depth);
            set(ax, 'YDir', 'normal')
            title("Before Chopping")
            axis tight equal
            ylabel("Depth [mm]")
            xlabel("Depth [mm]")
            ax = subplot(2,2,3);
            imagesc(ax, this.grid.int.axialVec, this.grid.int.depthVec,  log(abs(this.data.int.depth)));
            set(ax, 'YDir', 'normal')
            title("Before Chopping [Log]")
            axis tight equal
            ax = subplot(2,2,2);
            imagesc(ax, this.grid.spat.axialVec, this.grid.spat.depthVec,  this.data.spat.depth);
%             imagesc(ax, this.data.spat.depth);
            axis tight equal
            set(ax, 'YDir', 'normal')
            title("After Chopping")
            ax = subplot(2,2,4);
            imagesc(ax, this.grid.spat.axialVec, this.grid.spat.depthVec, log(abs(this.data.spat.depth)));
            axis tight equal
            set(ax, 'YDir', 'normal')
            title("After Chopping [Log]")


            figure();
            subplot(1,3,1)
            imagesc(this.grid.pulses.axialVec, this.grid.pulses.pulseStartPos, this.data.pulses.pulses)
            axis tight equal;
            xlabel("Depth [mm]")
            ylabel("Distance From Transducer [mm]")
            subplot(1,3,2)
            imagesc(this.grid.pulses.axialVec, this.grid.pulses.depthVec, this.data.pulses.pulsesEnv)
            axis tight equal;
            xlabel("Depth [mm]")
            ylabel("Distance From Transducer [mm]")
            subplot(2,3,3)
            plot(this.grid.pulses.axialVec, this.data.pulses.focalPulse)
            title("Pulse at Focal Point")
            xlabel("Depth[mm]")
            ylabel("Pressure [KPa]")
            subplot(2,3,6)
            plot(this.grid.pulses.axialVec, this.data.pulses.focalPulseEnv)
            title("Pulse Envelope at Focal Point")
            xlabel("Depth[mm]")
            ylabel("Pressure [KPa]")
        end
    end
end

