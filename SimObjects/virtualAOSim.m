 classdef virtualAOSim < handle
    %CREATEVIRTUALAOSIGNAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fullLoad
        
        usa;
        
        debug;
        displayDebug;
        debugTime;
        
        usVars;
        us;
        profiles;
        spatMat;
        pulses;
        vars;
        phiMath;
        
        figs;
    end
    
    methods (Static)
        function uVars = createUserVars()
            uVars.debug = false;
            uVars.displayDebug = false;
            uVars.debugTime = false;
            
            uVars.N   = 215;
            uVars.fUS = 1.25e6;
            uVars.pulseType = 'measured';
            
            uVars.spacerLen      = 6.4;
            uVars.spacerMaterial = 'PDMS';
            
            uVars.usFocalPointDistFromInterface = 30; %[mm]
            
            uVars.muEffVec = [0.074; 0.106; 0.152; 0.219; 0.323];
        end
    end
    
    methods
        function this = virtualAOSim()
            this.fullLoad     = false;
            
            this.debug        = false;
            this.displayDebug = false;
            this.debugTime    = false;
        end
        
        function loadFullUSPulse(this)
            fprintf("VAOS: Loading and Analyzing US Field - take a coffee break\n");
            this.usa = usAnalysis();

            uVars.usDataPath = 'AO_Transducer.mat';
            uVars.usDataType = '3D';
            uVars.intFactor = [1,1,1];
            this.usa.setVars(uVars);
            res = this.usa.analyse();
            usVars = this.usa.getVars();
            this.usVars.raw = usVars; %#ok<*PROP>
            % usa.cutPulses()
            % usa.calcPulsesProfile();
            % usa.extractFocalSignal();

            resSave.focalProfile      = res.focalProfile;
            resSave.focalPulsesEnvCut = res.focalPulsesEnvCut;
            resSave.focalSig          = res.focalSig;
            resSave.focalSigEnv       = res.focalSigEnv;
            resSave.focalPulseRawAx   = this.usVars.depthVec;
            
            this.us = resSave;
            save("C:\Users\sahiadl.EED\OneDrive - Technion\Graduate\Simulations\AcoustoOpticMuEffSimulation\analysedFocusedUS.mat", 'resSave', 'usVars', '-v7.3');
        end
        
        function loadUSPulse(this)
            fprintf("VAOS: Loading US Field\n");
            if this.fullLoad
                this.loadFullUSPulse();
            else
                res = load("C:\Users\sahiadl.EED\OneDrive - Technion\Graduate\Simulations\AcoustoOpticMuEffSimulation\analysedFocusedUS.mat");
                this.usVars.raw = res.usVars;
                this.us         = res.resSave;
                this.us.focalProfileMat = this.us.focalProfile;
            end            
        end

        function setVars(this, uVars)
            fprintf("VAOS: Set User Vars\n");
            this.debug = uVars.debug;
            this.displayDebug = uVars.displayDebug;
            this.debugTime = uVars.debugTime;
            
            this.vars.N =  uVars.N;
            this.vars.fUS = uVars.fUS;
            this.vars.pulseType = uVars.pulseType;
            this.vars.spacerLen = uVars.spacerLen;
            this.vars.spacerMaterial = uVars.spacerMaterial;
            this.vars.usFocalPointDistFromInterface = uVars.usFocalPointDistFromInterface; %[mm]
            
            this.vars.muEffVec = uVars.muEffVec;
            this.vars.numMu = length(this.vars.muEffVec);
        end
        
        function simVars = getVars(this)
           simVars = this.vars; 
        end
        
        function extractUSVars(this)
            fprintf("VAOS: Extracting US Profile and Variables\n");
            this.us.focalProfile     = squeeze(this.us.focalProfileMat(9,9,:));
            this.us.focalProfileNorm = this.us.focalProfile/max(this.us.focalProfile);
            this.us.depthProfile     = envelope(this.us.focalProfile, 25, 'peaks');
            [~, this.usVars.profileFocalPointIdx] = max(this.us.depthProfile);
            this.us.depthProfileNorm = this.us.depthProfile/max(this.us.depthProfile);

            this.us.pulse = squeeze(this.us.focalPulsesEnvCut(9,9, this.usVars.raw.focalIndIntAx,:))';
            this.us.pulseNorm = normMatf(this.us.pulse);
            this.usVars.pulseAx = this.usVars.raw.pulseVec;

            this.us.focalPulseSig = flip(this.us.focalSig(649:649+149-1));
            this.us.focalPulseSigNorm = this.us.focalPulseSig/max(this.us.focalPulseSig);

            this.usVars.profileDepthVecRaw = this.usVars.raw.pulseStartPosVec;
            this.usVars.usFocalLen = this.usVars.profileDepthVecRaw(this.usVars.profileFocalPointIdx);
            
            this.usVars.c  = round(this.usVars.raw.c);
            this.usVars.dX = this.usVars.raw.dAx;
            
            this.vars.tVec = this.usVars.raw.tVec;
            this.vars.dt   = this.vars.tVec(2) - this.vars.tVec(1);
        end
        
        function useCustomUSParams(this, use)
            if use
                this.vars.us.c = 1600;
                this.vars.dX   = this.vars.dt*this.vars.us.c*1e3;
                this.vars.us.profileDepthVecRaw = this.usVars.profileDepthVecRaw(1) + this.vars.dX*(0:1:length(this.usVars.profileDepthVecRaw) ) ;
                this.vars.us.usFocalLen = this.vars.us.profileDepthVecRaw(this.usVars.profileFocalPointIdx);
            else
                this.vars.us.c = this.usVars.c;
                this.vars.dX   = this.usVars.dX;
                this.vars.us.profileDepthVecRaw = this.usVars.profileDepthVecRaw;
                this.vars.us.usFocalLen         = this.usVars.usFocalLen;
            end
        end
        
        function calcSimDimension(this)
            fprintf("VAOS: Calculating Simulation Dimensions\n");
            this.vars.singleCycleLen    = (1/this.vars.fUS)*this.vars.us.c*1e3;
            this.vars.singleCycleIdxRaw = this.vars.singleCycleLen/this.vars.dX;
            this.vars.singleCycleIdx    = round(this.vars.singleCycleIdxRaw);
            
            if this.debug
                this.vars.sqncCompSpac = 8;
            else
                this.vars.sqncCompSpac = this.vars.singleCycleIdx;
            end
            
            this.vars.reconSize = this.vars.N*this.vars.sqncCompSpac;
        end
        
        function createSpace(this)
            fprintf("VAOS: Creating Spatial Vars:\n");
            dX = this.vars.dX;
            reconSize = this.vars.reconSize;
            
            % Creating system axis:
            xRaw = (0:1:reconSize-1)*dX;
            x = xRaw - (this.vars.us.usFocalLen + this.vars.usFocalPointDistFromInterface); % mm
            x = xRaw - mean(xRaw); % mm
            [~, minIdx] = min(abs(x));
            x = x - x(minIdx); % to ensure '0' exist on the grid.
            
            if this.debug
                % this is to make sure that even with short debug
                % simulation we still simulate around the '0' point.
                x = x - x(round(length(x)/2));
            end
            
            [xOrigin, xOriginIdx] = min(abs(x));

            % Incident:           
            incIdxStart = 1;
            incIdxEnd   = xOriginIdx;
            incIdx      = incIdxStart:incIdxEnd;
            
            incLenIdx = length(incIdx);

            incXStart = x(incIdxStart);
            incXEnd   = x(incIdxEnd); %  0
            incX      = x(incIdxStart:incIdxEnd);
            
            % Spacer:
            spacerLenEff = this.vars.spacerLen*2;

            spacerX      = xOrigin + (0:dX:spacerLenEff) + dX;
            spacerXStart = spacerX(1);
            spacerXEnd   = spacerX(end);
            
            spacerLenIdx = length(spacerX);
            
            spacerIdx      = xOriginIdx + (1:1:spacerLenIdx);
            spacerIdxStart = spacerIdx(1);
            spacerIdxEnd   = spacerIdx(end);
            
            spacerAirInt   = floor(spacerLenIdx/2) + spacerIdxStart;
            spacerIncIdx = spacerIdxStart:(spacerAirInt-1);
            spacerRedIdx = spacerAirInt:spacerIdxEnd;
            
            
            % Reflection space:
            refIdxStart = spacerIdxEnd+1;
            refIdxEnd   = reconSize;
            refIdx      = refIdxStart:refIdxEnd;
            
            refLenIdx = length(refIdx);
            
            refXStart = x(refIdxStart);
            refXEnd   = x(refIdxEnd);
            refX      = x(refIdxStart : refIdxEnd);
            
            % Collect parameters:
            this.vars.x = x;
            this.vars.space.xOriginIdx = xOriginIdx;
            
            this.vars.space.incIdxStart = incIdxStart;
            this.vars.space.incIdxEnd   = incIdxEnd;
            this.vars.space.incIdx      = incIdx;          
            this.vars.space.incLenIdx   = incLenIdx;
            this.vars.space.incXStart   = incXStart;
            this.vars.space.incXEnd     = incXEnd; %  0
            this.vars.space.incX        = incX;
            
            this.vars.space.spacerX        = spacerX;
            this.vars.space.spacerXStart   = spacerXStart;
            this.vars.space.spacerXEnd     = spacerXEnd;
            this.vars.space.spacerIdx      = spacerIdx;
            this.vars.space.spacerLenIdx   = spacerLenIdx;
            this.vars.space.spacerIdxStart = spacerIdxStart;
            this.vars.space.spacerIdxEnd   = spacerIdxEnd;
            this.vars.space.spacerAirInt   = spacerAirInt;
            this.vars.space.spacerIncIdx   = spacerIncIdx;
            this.vars.space.spacerRedIdx   = spacerRedIdx;
            
            this.vars.space.refIdxStart = refIdxStart;
            this.vars.space.refIdxEnd   = refIdxEnd;
            this.vars.space.refIdx      = refIdx;
            this.vars.space.refLenIdx   = refLenIdx;
            this.vars.space.refXStart   = refXStart;
            this.vars.space.refXEnd     = refXEnd;
            this.vars.space.refX        = refX;
        end
            
        function createMathematicalFluence(this)
            fprintf("VAOS: Creating Mathematical Fluence Profile\n");
            incX = this.vars.space.incX;
            refX = this.vars.space.refX;
            refXStart = this.vars.space.refXStart;
            spacerLenIdx = this.vars.space.spacerLenIdx;
            
            muEffVec = this.vars.muEffVec;
            numMu    = this.vars.numMu;

            phiInc = exp(-muEffVec.*abs(incX));
            spacer = zeros(numMu, spacerLenIdx);
            phiRef = 0*exp(-muEffVec.*abs(refX-refXStart));
            
            phi = [phiInc, spacer, phiRef];

            if this.debug
                reconSize = this.vars.reconSize;
                phi = (1e-2)*ones(numMu, reconSize);
                phi(:,1) = 1e-3;
                phi(:, 300:330) = 1;
            end
            
            this.phiMath = phi;
        end
        
        function alignAndInterpUS(this)
            fprintf("VAOS: Aligning US profile to Fluence\n");
            x = this.vars.x;
            
            [~, focalPointIdx] = min(abs(x+this.vars.usFocalPointDistFromInterface));
            usFocalPos = x(focalPointIdx);
            
            profileDepthVecArt = this.usVars.profileDepthVecRaw - this.usVars.usFocalLen + usFocalPos;
            depthProfileArt    = this.us.depthProfile';
            
            if x(1) < profileDepthVecArt(1)
                profileDepthVecArt = [x(1), profileDepthVecArt];
                depthProfileArt = [0, this.us.depthProfile'];
            end
            
            if x(end) > profileDepthVecArt(end)
                profileDepthVecArt = [profileDepthVecArt, x(end)];
                depthProfileArt    = [depthProfileArt, 0];
            end
            
            focalProfileInt = normMatf(interp1(profileDepthVecArt, depthProfileArt, x, 'pchip'));
            
            this.vars.focalPointIdx = focalPointIdx;
            this.profiles.focalProfileInt = focalProfileInt;
        end
        
        function createAcousticInterface(this)
            fprintf("VAOS: Creating Acoustic Reflection\n");
            reconSize = this.vars.reconSize;
            incIdx = this.vars.space.incIdx;
            spacerIncIdx = this.vars.space.spacerIncIdx;
            spacerRedIdx = this.vars.space.spacerRedIdx;
            refIdx = this.vars.space.refIdx;

            
            zPDMS  = 1.048e6;   %[Ns * m^-3]
            zWater = 1.494e6;   %[Ns * m^-3]
            zAir   = 1.2 * 343; % [Ns * m^-3] 
            
            switch this.vars.spacerMaterial
                case 'PDMS'
                    R12 = abs((zPDMS - zWater) / (zPDMS + zWater));
                    T12 = abs( 2*zPDMS / (zPDMS + zWater));
                    R23 = abs((zAir - zPDMS) / (zAir + zPDMS));
                    T21 = abs( 2*zWater / (zPDMS + zWater));

                case 'Water'
                    R12 = 0;
                    T12 = 1;
                    R23 = abs((zAir - zWater) / (zAir + zWater));
                    T21 = 1;
            end
            
            usIntProfile = zeros(1, reconSize); 
            usIntProfile (incIdx)       = 1;
            usIntProfile (spacerIncIdx) = T12;
            usIntProfile (spacerRedIdx) = T12*R23;
            usIntProfile (refIdx)       = T12*R23*T21;

            this.profiles.usIntProfile = usIntProfile;
        end
        
        function naiveSim(this)
            fprintf("VAOS: Performing Naive Simulation\n");
            measPhiNaive = this.phiMath.*this.profiles.focalProfileInt.*this.profiles.usIntProfile;

            figure();
            for i=1:5
                ax(i) = subplot(2,3,i); %#ok<AGROW>
                hold on
                plot(log(this.phiMath(i,:)))
                plot(real(log(measPhiNaive(i,:))))
                xlabel("X[mm]")
                ylabel("Fluence [AU]")
                title(sprintf("Math. Phantom: %d", i))
            end
            linkaxes(ax);
        end
        
        function createPulses(this)
            fprintf("VAOS: Creating Pulses\n");
            singleCycleIdx = this.vars.singleCycleIdx;
            sqncCompSpac = this.vars.sqncCompSpac;
            
            % Ideal Delta:
            this.pulses.pulseIdeal    = [zeros(1,148), 1];
            this.pulses.pulseIdealSig = this.pulses.pulseIdeal;
            
            % Window:
            this.pulses.pulseWindow    = [zeros(1,149-singleCycleIdx), ones(1,singleCycleIdx)];
            this.pulses.pulseWindowSig = this.pulses.pulseWindow;
            
            % Triangle:
            this.pulses.pulseTriangle    = [zeros(1,149-singleCycleIdx), linspace(0,1,singleCycleIdx)];
            this.pulses.pulseTriangleSig = this.pulses.pulseTriangle;
            
            % Ideal Single Cycle Sinus:
            pad = sqncCompSpac-singleCycleIdx;
            
            sinSig = sin(pi*((1:1:(singleCycleIdx))/singleCycleIdx)).^2;
            this.pulses.pulseSin = [zeros(1,pad), sinSig];
            
            sinSig = -sin(2*pi*((1:1:(singleCycleIdx))/singleCycleIdx));
            this.pulses.pulseSinSig = [zeros(1,pad), sinSig];

            % Measured Pulse:
            this.pulses.pulseMeas    = this.us.pulseNorm;
            this.pulses.pulseMeasSig = this.us.focalPulseSigNorm;
        end
        
        function [env, sig] = choosePulseShape(this)
            switch this.vars.pulseType
                case 'delta'
                    env = this.pulses.pulseIdeal;
                    sig = this.pulses.pulseIdealSig;
                case 'sine'
                    env = this.pulses.pulseSin;
                    sig = this.pulses.pulseSinSig;
                case 'window'
                    env = this.pulses.pulseWindow;
                    sig = this.pulses.pulseWindowSig;
                case 'triangle'
                    env = this.pulses.pulseTriangle;
                    sig = this.pulses.pulseTriangleSig;
                case 'measured'
                    env = this.pulses.pulseMeas;
                    sig = this.pulses.pulseMeasSig;
                case 'other'   
            end
        end
        
        function buildSPMatrix(this)
            fprintf("VAOS: Creating SP Matrix\n");
            reconSize       = this.vars.reconSize;
            focalProfileInt = this.profiles.focalProfileInt;
            usIntProfile    = this.profiles.usIntProfile;
            if this.debug
                focalProfileInt = ones(1, reconSize);
                usIntProfile    = ones(1, reconSize);
            end
            
            [pulseEnv, pulseSig] = this.choosePulseShape();
            pulseLen = length(pulseEnv);
            
            padSize  = pulseLen-1;
            endIdx   = reconSize;
            startIdx = 1;

            if this.displayDebug
                hFig1 = figure();
                subplot(3,3,1); h1 = plot(zeros(1,reconSize));   title("CurSqnc");           ylim([-1,1]);
                subplot(3,3,2); h2 = plot(zeros(1,2*reconSize)); title("curSqncSpatial2");   ylim([-1,1]);
                subplot(2,3,3); h3 = imagesc(zeros(reconSize));  title("transMatEnvSP");     axis tight equal; colorbar;
                subplot(3,3,4); h4 = plot(zeros(1,2*reconSize)); title("transmissionEnvSP"); ylim([-1,1]);
                subplot(3,3,5); h5 = plot(zeros(1,reconSize));   title("transMatEnvSP");     ylim([-1,1]);
                subplot(2,3,6); h6 = imagesc(zeros(reconSize));  title("transMatEnvSP");     axis tight equal; colorbar;
                subplot(3,3,7); h7 = plot(zeros(1,2*reconSize)); title("transmissionSigSP"); ylim([-1,1]);
                subplot(3,3,8); h8 = plot(zeros(1,reconSize));   title("transMatSigSP");     ylim([-1,1]);
            end

            transMatEnvSP = zeros(reconSize);
            transMatSigSP = zeros(reconSize);
            
            pulsePos    = [1, zeros(1,reconSize-1)];
            for i = 1:reconSize
                curSqnc          = circshift(pulsePos, (i-1));
                curSqncSpatial   = curSqnc.*focalProfileInt.*usIntProfile;
                transmissionEnvSP  = conv(curSqncSpatial, flip(pulseEnv), 'full');
                transMatEnvSP(:,i) = transmissionEnvSP(startIdx:endIdx);

                transmissionSigSP  = conv(curSqncSpatial, flip(pulseSig), 'full');
                transMatSigSP(:,i) = transmissionSigSP(startIdx:endIdx);

                if this.displayDebug && ((i==1) || (~mod(i,10)) || (i ==reconSize))
                    set(h1, 'YData', curSqnc);
                    set(h2, 'YData', curSqncSpatial);
                    set(h3, 'CData', transMatEnvSP);
                    set(h4, 'YData', transmissionEnvSP);
                    set(h5, 'YData', transMatEnvSP(i,:));
                    set(h6, 'CData', transMatSigSP);
                    set(h7, 'YData', transmissionSigSP);
                    set(h8, 'YData', transMatSigSP(i,:));
                    drawnow();
                end
            end
            this.spatMat.transMatEnvSP     = transMatEnvSP;
            this.spatMat.transMatSigSP     = transMatSigSP;
            this.spatMat.transMatEnvNormSP = transMatEnvSP/max(transMatEnvSP(:));
            this.spatMat.transMatSigNormSP = transMatSigSP/max(transMatSigSP(:));
        end
        
        function res = reconSP(this, phi, figs)
            % Envelope Reconstruction
            phiEnvRecon     = this.spatMat.transMatEnvNormSP * phi';
            phiEnvReconNorm = phiEnvRecon ./ max(phiEnvRecon,[],1);
            
            % Signal Reconstruction (meaningless)
            phiSigRecon     = this.spatMat.transMatSigNormSP * phi';
            phiSigReconNorm = phiSigRecon ./ max(phiSigRecon,[],1);

            res.phiEnvRecon = phiEnvRecon;
            res.phiEnvReconNorm = phiEnvReconNorm;
            res.phiSigRecon = phiSigRecon;
            res.phiSigReconNorm = phiSigReconNorm;
            
            if figs
                hFig1 = figure();
                set(hFig1, 'NumberTitle', 'off', 'Name', 'SP Recon Envelope');
                for i=1:5
                    ax(i) = subplot(2,3,i); %#ok<AGROW>
                    hold on
                    plot(this.vars.x, log(phi(i,:)))
                    plot(this.vars.x, real(log(phiEnvReconNorm(:,i))))
                    xlabel("X[mm]")
                    ylabel("Fluence [AU]")
                    title(sprintf("Phantom: %d", i))
                    xlim([-60,30]);
                    ylim([-10,0]);
                end
                linkaxes(ax);
                
%                 hFig2 = figure();
%                 set(hFig2, 'NumberTitle', 'off', 'Name', 'SP Recon Sig');
%                 for i=1:this.vars.numMu
%                     ax(i) = subplot(2,3,i);
%                     hold on
%                     plot( log(phi(i,:)))
%                     plot( real(log(phiSigReconNorm(:,i))))
%                     xlabel("X[mm]")
%                     ylabel("Fluence [AU]")
%                     title(sprintf("Phantom: %d", i))
%                     ylim([-10,0])
%                     xlim([550, 1200])
%                 end
%                 linkaxes(ax);
            end
        end
        
        function buildHadMatrix(this)
            fprintf("VAOS: Creating Had Matrix\n");
            N = this.vars.N;
            reconSize       = this.vars.reconSize;
            sqncCompSpac    = this.vars.sqncCompSpac;
            focalProfileInt = this.profiles.focalProfileInt;
            usIntProfile    = this.profiles.usIntProfile;
            if this.debug
                focalProfileInt = ones(1, reconSize);
                usIntProfile    = ones(1, reconSize);
            end
            
            % Retreive Pulse shape:
            [pulseEnv, pulseSig] = this.choosePulseShape();
            
            % Create Matrix Construction:
            this.vars.sMat   = createSMatrix(N);
            sVec   = this.vars.sMat(1, :);
            this.vars.sVec = sVec;
            
            sVecInt          = zeros(sqncCompSpac,N);
            sVecInt(1, :)  = this.vars.sVec;
            sVecInt          = sVecInt(:)';
            
            % Calculate Convolution indices:
            padSize = length(pulseEnv)-1;
            transFullLen = 3*reconSize+padSize;
            startIdx = reconSize + 1;
            endIdx   = 2*reconSize;
%             startIdx = reconSize + padSize + 1;
%             endIdx   = 2*reconSize + padSize;
            
            idxVec = 1:transFullLen;
            frame  = (idxVec >= startIdx) & (idxVec <= endIdx);
            
            if this.displayDebug
                figure()
                ax(1) = subplot(3,2,1); h1 = plot(zeros(1,reconSize));    title("CurSqnc");
                ax(3) = subplot(2,2,2); h3 = imagesc(zeros(reconSize));   title("transMatEnvHad"); axis tight equal; colorbar;
                ax(5) = subplot(3,2,3); h5 = plot(zeros(1,reconSize));    title("transMat");
                ax(6) = subplot(2,2,4); h6 = imagesc(zeros(reconSize));   title("transMatEnvHad"); axis tight equal; colorbar;
                ax(8) = subplot(3,2,5); h8 = plot(zeros(1,reconSize));    title("transMat");
            end

            transEnvMat = zeros(reconSize);
            transSigMat = zeros(reconSize);

            idxVec = 1:reconSize;
            
            transMatEnvSP = this.spatMat.transMatEnvNormSP;
            transMatSigSP = this.spatMat.transMatSigNormSP;
            
            for i = 1:reconSize               
                curSqnc = circshift(sVecInt, (i-1));
                curIdxs = idxVec(logical(curSqnc));
                transEnvMat(i,:) = sum(transMatEnvSP(curIdxs,:),1);
                transSigMat(i,:) = sum(transMatSigSP(curIdxs,:),1);

                if this.displayDebug  && ((i==1) || (~mod(i,10)) || (i ==reconSize))
                    set(h1, 'YData', curSqnc);
                    set(h3, 'CData', transEnvMat);
                    set(h5, 'YData', transEnvMat(i,:));
                    set(h6, 'CData', transSigMat);
                    set(h8, 'YData', transSigMat(i,:));
                    drawnow();
                end
            end
            
            this.spatMat.transMatEnvHad     = transEnvMat;
            this.spatMat.transMatSigHad     = transSigMat;
            this.spatMat.transMatEnvNormHad = transEnvMat/ max(transEnvMat(:));
            this.spatMat.transMatSigNormHad = transSigMat/ max(transSigMat(:));
        end
        
        function buildHadInvMat(this)
            N            = this.vars.N;
            sqncCompSpac = this.vars.sqncCompSpac;
            reconSize    = this.vars.reconSize;
            
            this.vars.sMatInv = inv(this.vars.sMat);
            sVecInv = this.vars.sMatInv(1,:);

            sVecInvSqnc        = zeros(sqncCompSpac,N);
            sVecInvSqnc(1,:)   = flip(sVecInv);
            sVecInvSqnc        = sVecInvSqnc(:)';

            sMatInvSqnc = zeros(reconSize);

            if this.displayDebug
                figure()
                h1 = imagesc(zeros(reconSize)); title("Inverse Hadamard Matrix");
            end
            
            for i=1:reconSize
                sMatInvSqnc(i,:) = circshift(sVecInvSqnc, -(i-1));
                if this.displayDebug && ~mod(i,20)
                   set(h1, 'CData',  sMatInvSqnc);
                   drawnow();
                end
            end

            this.spatMat.sMatInvSqnc = flip(sMatInvSqnc,1);
        end
        
        function res = reconHad(this, phi, figs)
            phiEnvRecon     = this.spatMat.sMatInvSqnc * this.spatMat.transMatEnvNormHad * phi';
            phiEnvReconNorm = phiEnvRecon./max(phiEnvRecon,[],1);
            
            phiSigRecon     = this.spatMat.sMatInvSqnc * (this.spatMat.transMatSigNormHad * phi');
            phiSigReconNorm = phiSigRecon./max(phiSigRecon,[],1);
            
            res.phiEnvRecon     = phiEnvRecon;
            res.phiEnvReconNorm = phiEnvReconNorm;
            res.phiSigRecon     = phiSigRecon;
            res.phiSigReconNorm = phiSigReconNorm;
            
            if figs
                hFig1 = figure();
                set(hFig1, 'NumberTitle', 'off', 'Name', 'Had Recon Sig');
                for i=1:this.vars.numMu
                    ax(i) = subplot(2,3,i);
                    hold on
                    plot(this.vars.x, log(phi(i,:)))
                    plot(this.vars.x, log(abs(phiEnvReconNorm(:,i))))
%                     plot(log(phi(i,:)))
%                     plot(log(abs(phiEnvReconNorm(:,i))))
%                     plot(phi(i,:))
%                     plot(phiEnvReconNorm(:,i))
                    xlabel("X[mm]")
                    ylabel("Fluence [AU]")
                    title(sprintf("Phantom: %d", i))
                    xlim([-60,30]);
                    ylim([-10,0]);
                end   
                linkaxes(ax);
%                 hFig2 = figure();
%                 set(hFig2, 'NumberTitle', 'off', 'Name', 'Had Recon Sig');
%                 for i=1:this.vars.numMu
%                     subplot(2,3,i)
%                     hold on
%                     plot(this.vars.x, log(phi(i,:)))
%                     plot(this.vars.x, log(abs(phiSigReconNorm(:,i))))
%                     xlabel("X[mm]")
%                     ylabel("Fluence [AU]")
%                     title(sprintf("Phantom: %d", i))
%     %                 xlim([-60,30]);
%     %                 ylim([-10,0]);
%                 end  
            end
            
        end
        
        function res = reconAll(this, phi)
            res.sp = this.reconSP(phi, false);
            res.had = this.reconHad(phi, false);
            
            hFig1 = figure();
            set(hFig1, 'NumberTitle', 'off', 'Name', 'Had Recon Sig');
            for i=1:this.vars.numMu
                ax(i) = subplot(2,3,i);
                hold on
                plot(this.vars.x, log(phi(i,:)))
                plot(this.vars.x, log(abs(res.sp.phiEnvReconNorm(:,i))))
                plot(this.vars.x, log(abs(res.had.phiEnvReconNorm(:,i))))
                xlabel("X[mm]")
                ylabel("Fluence [AU]")
                title(sprintf("Phantom: %d", i))
                legend("Phi", "SP", "Had", 'Location', 'northwest')
                xlim([-60,30]);
                ylim([-10,0]);
            end   
            linkaxes(ax);
            
        end
        
        function phiSimAligned = alignExternalPhi(this, phi, depthVec, figs)
            dX = this.vars.dX;
            x = this.vars.x;
            numMu = this.vars.numMu;
            spacerLen = this.vars.space.spacerLenIdx;

            % Interpolate for x vec resolution:
            depthVecSimInt = depthVec(1) : dX : depthVec(end);
            phiSim1SideInt = zeros(numMu, length(depthVecSimInt));
            phiSim1SideLen = length(depthVecSimInt);
            
            for i=1:numMu
                phiSim1SideInt(i,:) = interp1(depthVec, phi(i,:), depthVecSimInt, 'pchip');
            end
            
            % Build spatial fluence:
            phiSim2Side   = [flip(phiSim1SideInt,2), zeros(numMu,spacerLen), phiSim1SideInt];
            depthVec2Side = (0:1:(size(phiSim2Side,2)-1))*dX;
            xSim          = depthVec2Side - depthVec2Side(phiSim1SideLen);
            
            if x(1) < xSim(1)
                phiSim2Side = [zeros(numMu,1), phiSim2Side];
                xSim = [x(1), xSim];
            end
            
            if x(end) > xSim(end)
                phiSim2Side    = [phiSim2Side, zeros(numMu,1)];
                xSim = [xSim, x(end)];
            end
            
            
            phiSimAligned = zeros(5, length(x));
            for i=1:5
                curSim = interp1(xSim, phiSim2Side(i,:), x, 'pchip');
                phiSimAligned(i,:) = curSim/max(curSim);
            end
            
            if figs
                figure();
                subplot(1,2,1)
                plot(x, phiSimAligned)
                subplot(1,2,2)
                plot(x, log(normMatf(phiSimAligned,2)))
            end
        end
        
        function matchMeasAndSim(this, phiMeas, xMeas, phiSim, xSim)
            phiRawAlign = this.alignExternalPhi(phiSim, xSim, false);
            phiReconHad = this.reconHad(phiRawAlign, false);
            phiReconHadNorm = phiReconHad.phiEnvReconNorm;
            
            figure();
            for i=1:this.vars.numMu
                ax(i) = subplot(2,3,i);
                hold on
                plot(this.vars.x, log(phiRawAlign(i,:)))
                plot(this.vars.x, log(abs(phiReconHadNorm(:,i))))
                plot(xMeas, log(phiMeas(i,:)));
                xlabel("X[mm]")
                ylabel("Fluence [AU]")
                legend("MCX", "Conv", "Meas", 'Location', 'northwest')
                title(sprintf("Phantom: %d", i))
                xlim([-50,50])
                ylim([-10,0])
            end
            linkaxes(ax)
            
        end

        function res = speckleSim(this, phi, figs)
            phi = gpuArray(abs(phi));
            numPh = size(phi,1);

            %----------------------------------------
            % Calculate AO Reconstruction dimensions:
            %----------------------------------------
            reconSize       = this.vars.reconSize;
            
            framesPerSig    = 1000; % number of speckle-decorrelation frames
            sqncPerFrame    = 10;
            pulsePerSqnc    = this.vars.N;
            samplesPerPulse = this.vars.singleCycleIdx; % theoretical
            
            samplesPerSqnc  = samplesPerPulse * pulsePerSqnc;
            samplesPerFrame = samplesPerSqnc  * sqncPerFrame;
            samplesPerPos   = samplesPerPulse * sqncPerFrame;
            
            numOfPos = pulsePerSqnc;
            
            
            batchSize    = 10;
            numOfBatch   = ceil(framesPerSig/batchSize);
            framesPerSig = numOfBatch * batchSize;
            %----------------------------------------
            % Calculate Time-Space Grid:
            %----------------------------------------
            x = this.vars.x;
            xMat = repmat(gpuArray(x), reconSize, 1, batchSize); % Time x Space
%             xMat = repmat(gpuArray(x), reconSize, 1); % Time x Space
            xUS = x(1:samplesPerPulse:end);
            this.vars.xUS = xUS;
            
            dt  = this.vars.dt; 
            fs  = 1/dt;
            fUS = this.vars.fUS;
            
            N = samplesPerPos; % only for readability of the following code:
            fBar      = (fs/N) *  ( (-N/2) : 1 : (N/2)-1 ) *1e-6;
            df        = fs / N;
            fUsIdxPos = floor(fUS /df) + 1 + N/2;
            
            %----------------------------------------
            % EM Parameters:
            %----------------------------------------
            transMatSigNormSP  = this.spatMat.transMatSigNormSP; % US Modulation
            transMatSigNormHad = this.spatMat.transMatSigNormHad;
            
            n      = 1.34;
            lambda = 785e-9;
            k0     = 2*pi/lambda; 
            gamma  = 1e-1; % modulation depth
            dnSP   = gpuArray(-pi*gamma*transMatSigNormSP); % IS dn linear? should find a reference
            dnHad  = gpuArray(-pi*gamma*transMatSigNormHad);
            
            numOfGrain = 100;  % number of speckle grains;
            SBR        = 1000; % Signal-to-background ratio -> does not affect the SNR (obviously)
            
            cleanSigSP = gpuArray(zeros(samplesPerSqnc, framesPerSig));
            IdSPMat  = zeros(samplesPerFrame, framesPerSig, numPh);
            IdHadMat = zeros(samplesPerFrame, framesPerSig, numPh);
            
            %-----------------------------------
            % Create Electrical Circuit LPF:
            %-----------------------------------
            [~,  lpf]  = lowpass(zeros(1,reconSize), 2*fUS, fs, 'StopbandAttenuation', 40);
            lpfT       = lpf.Coefficients;
            lpfLen     = length(lpfT);
            ipfPivot   = floor(lpfLen/2)+1;
            padSize    = ipfPivot-1;
            paddedSize = samplesPerFrame + 2*padSize;
            
            filterMat  = zeros(samplesPerFrame, paddedSize, 'gpuArray');
            
            % Create a batch:
            batchRowSize    = 1000;
            batchColSize = batchRowSize + 2*padSize;
            batchMat     = zeros(batchRowSize, batchColSize, 'gpuArray');
            
            for i=1:batchRowSize
               endCol   = i + lpfLen-1;
               batchMat(i,i:endCol) = lpfT;
            end
            
            numBatch    = floor(samplesPerFrame/batchRowSize);
            batchedRows = batchRowSize*numBatch;
            
            % Insert batch to filter matrix:
            for i = 1:numBatch
                startRowBatch = (i-1)*batchRowSize+1;
                endRowBatch   = i*batchRowSize;
                startColBatch = (i-1)*batchRowSize+1;
                endColBatch   = (i-1)*batchRowSize + batchColSize;
                filterMat(startRowBatch:endRowBatch, startColBatch:endColBatch) = batchMat;
            end
            
            % Complete the residual lines:
            for i = batchedRows+1 : samplesPerFrame
                endCol   = i + lpfLen-1;
                filterMat(i,i:endCol) = lpfT;
            end
            
            clear batchMat
            %-----------------------------------
            % Create Inverse Hadamard Matrix:
            %-----------------------------------
%             sVecInv = this.vars.sMatInv(1,:);
% 
%             sVecInvSqnc        = zeros(sqncCompSpac,N);
%             sVecInvSqnc(1,:)   = flip(sVecInv);
%             sVecInvSqnc        = sVecInvSqnc(:)';
% 
%             sMatInvSqnc = zeros(reconSize);
%             
%             for i=1:reconSize
%                 sMatInvSqnc(i,:) = circshift(sVecInvSqnc, -(i-1));
%             end
%             
%             figure();
%             subplot(1,2,1)
%             imagesc(transMatSigNormHad)
%             subplot(1,2,2)
%             imagesc(sMatInvSqnc)
            
            sMatInvSqnc = this.spatMat.sMatInvSqnc;
%             sMatInv = eye(sqncPerFrame);
            
%             sMatInv = kron(eye(sqncPerFrame), sMatInvSqnc);
            sMatInv = gpuArray(sMatInvSqnc);
            %----------------------------------------
            
            phiReconSP      = zeros(numPh, numOfPos);
            phiReconHad     = zeros(numPh, numOfPos);
            phiReconSPNorm  = zeros(numPh, numOfPos);
            phiReconHadNorm = zeros(numPh, numOfPos);
            
            if this.displayDebug
                figure();
                subplot(1,2,1)
                plot(lpfT)
                subplot(1,2,2)
                imagesc(filterMat)
                
                figure()
                subplot(2,2,1)
                hp1 = plot(x, zeros(1,reconSize));
                title("Clean Sig")
                subplot(2,2,2)
                hp2 = plot(x, zeros(1,reconSize));
                title("Noised Signal")
                subplot(2,2,3)
                hp3 = plot(x, zeros(1,reconSize));
                title("Digitized Signal")
            end
            
            for j=1:numPh
                fprintf("SpeckleSim: Phantom #%d \n", j);
                Ephii = repmat(sqrt(phi(j,:)), reconSize, 1);
                Ebkg  = SBR*sum(sqrt(phi(j,:)))* ones(reconSize, numOfGrain);
                
                %----------------------------------------------------------
                % Intereferece between modulated and unmodulated fields:
                %----------------------------------------------------------
                trun = tic;
                tFrames = tic;
                for i=1:numOfBatch
                    tFrame = tic;
                    startIdx =  (i-1)*batchSize+1;
                    endIdx   =  i*batchSize;
                    
                    tPhase = tic;
                    phaseX = 2*pi*rand(1,numOfPos, batchSize, 'gpuArray');
                    phaseX = repmat (phaseX, samplesPerPulse, 1, 1);
                    phaseX = reshape(phaseX, 1, samplesPerSqnc, batchSize);
                    phaseX = repmat (phaseX, reconSize, 1, 1);
                    T(i, 1) = toc(tPhase);
                    
                    tField = tic;
                    EiSP  = Ephii .* exp(1i*n*k0.*xMat) .* exp(1i*phaseX) .* exp(1i*dnSP);
                    EiHad = Ephii .* exp(1i*n*k0.*xMat) .* exp(1i*phaseX) .* exp(1i*dnHad);
                    T(i, 2) = toc(tPhase);
                    
                    tGrains = tic;
                    EgAOSP  = repmat(sum(EiSP,2),1, numOfGrain, 1);
                    EgAOHad = repmat(sum(EiHad,2),1, numOfGrain, 1);
                    T(i, 3) = toc(tPhase);
                    
                    tBkg = tic;
                    phaseG  = repmat(2*pi*rand(1, numOfGrain, batchSize, 'gpuArray'), reconSize, 1, 1);
                    phaseG2 = exp(1i*phaseG);
                    Ebkgg   = Ebkg.*phaseG2;
                    T(i, 4) = toc(tPhase);
                    
                    tInt = tic;
                    EgSP  = EgAOSP  + Ebkgg;
                    EgHad = EgAOHad + Ebkgg;
                    T(i, 5) = toc(tPhase); 
                    
                    tPower = tic;
                    IgSP  = EgSP.*conj(EgSP);
                    IgHad = EgSP.*conj(EgHad);
                    T(i, 6) = toc(tPhase);
                    
                    % Calculate net signal that hits the detector:
                    tDet = tic;
                    cleanSigSP(:,startIdx:endIdx)  = squeeze(sum(abs(IgSP), 2));
                    cleanSigHad(:,startIdx:endIdx) = squeeze(sum(abs(IgHad), 2));
                    T(i, 7) = toc(tPhase);
                    
                    T(i, 8) = toc(tFrame);
                end
                avgT = 1e3*mean(T,1);
                fprintf("\nSpeckleSim: create random data for all frames %.2f[s]\n", toc(tFrames))
                fprintf ("SpeckleSim: Average Frame Calc. Timings:\n")
                fprintf ("Phase: %.2f | Field: %.2f | Grains: %.2f | Bkg: %.2f | Int: %.2f | Power: %.2f | Det: %.2f\n",...
                    avgT(1), avgT(2), avgT(3), avgT(4), avgT(5), avgT(6), avgT(7));
                
                tClearFrame = tic;
                clear Ephii Ebkg;
                clear phaseX EiSP EiHad EgAOSP EgAOHad phaseG phaseG2 Ebkgg EgSP EgHad IgSP IgHad 
                fprintf("SpeckleSim: Clearing Variables  %.2f[s]\n", toc(tClearFrame));
                
                % replicate the clean signal:
                tSqnc = tic;
                cleanSigFrameSP  = repmat(cleanSigSP,  sqncPerFrame, 1);
                cleanSigFrameHad = repmat(cleanSigHad, sqncPerFrame, 1);
                fprintf("SpeckleSim: Replicate clean signal for Sqnc %.2f[s]\n", toc(tSqnc))
                
%                 cleanSigFrameHadSqnc   = reshape(cleanSigFrameHad(:,1), samplesPerSqnc, sqncPerFrame);
%                 cleanSigFrameHadSqncDM = sMatInv*cleanSigFrameHadSqnc;
%                 cleanSigFrameHadDM     = reshape(cleanSigFrameHadSqncDM, samplesPerFrame, 1);
%                 
%                 figure();
%                 ax(1) = subplot(1,2,1);
%                 plot(cleanSigFrameSP(:,1));
%                 ax(2) = subplot(1,2,2);
%                 plot(cleanSigFrameHadDM(:,1));
%                 linkaxes(ax, 'x');
                
                %----------------------------------------------------------
                % Intereferece between modulated and unmodulated fields:
                %----------------------------------------------------------
                tNoise = tic;
                noise         = sqrt(mean(cleanSigSP(:))) * randn(samplesPerFrame, framesPerSig);
                noisedSigSP   = noise + cleanSigFrameSP;
                noisedSigHad  = noise + cleanSigFrameHad;
                fprintf("SpeckleSim: Add Noise %.2f[s]\n", toc(tNoise))
                
                % Signal after AC coupling:
                tAC = tic;
                acCoupledSP  = noisedSigSP  - mean(cleanSigFrameSP(:));
                acCoupledHad = noisedSigHad - mean(cleanSigFrameHad(:));
                fprintf("SpeckleSim: Ac Coupling %.2f[s]\n", toc(tAC));
                
                % Signal after LPF:
                tPad = tic;
                paddedSigSP  = [acCoupledSP(1)   * ones(padSize,framesPerSig);...
                                acCoupledSP; ...
                                acCoupledSP(end) * ones(padSize,framesPerSig)];
                paddedSigHad = [acCoupledHad(1)  * ones(padSize,framesPerSig);...
                                acCoupledHad; ...
                                acCoupledHad(end)* ones(padSize,framesPerSig)];
                fprintf("SpeckleSim: Padding %.2f[s]\n", toc(tPad));
                
                tFilter = tic;
                IdSP  = filterMat * paddedSigSP;
                IdHad = filterMat * paddedSigHad;
                fprintf("SpeckleSim: Filtering %.2f[s]\n", toc(tFilter))
                
                clear cleanSigFrameSP cleanSigFrameHad 
                clear noise noisedSigSP noisedSigHad acCoupledSP 
                clear acCoupledHad paddedSigSP paddedSigHad
                
                fprintf("SpeckleSim: Created the measured signal in %.2f[s]\n", toc(trun));
                
                %-----------------------------------
                % AOI Reconstruction Algorithm:
                %-----------------------------------
                tAOSP = tic;
                A0 = reshape(IdSP(:,:), samplesPerSqnc, sqncPerFrame, framesPerSig);
                A1 = permute(A0, [2,1,3]);
                A2 = reshape(A1, sqncPerFrame, samplesPerPulse, numOfPos, framesPerSig);
                A3 = permute(A2, [2,1,3,4]);
                A4 = reshape(A3, samplesPerPos, numOfPos, framesPerSig);                                   % AC-coupling for the entire signal
                A5 = A4 - repmat(mean(A4, 1), samplesPerPos, 1, 1);       % AC-Coupling for each channel - not sure needed when noise is present
                A6 = abs(fftshift(fft(A5,[],1),1)).^2;
                A7 = mean(A6, 3);
                curReconSP = gather(sqrt(A7(fUsIdxPos , :)));
                fprintf("SpeckleSim: ReconSP %.2f[s]\n", toc(tFilter));
                
                phiReconSP(j,:)  = curReconSP;
                
                maxVal = max(curReconSP);
                minVal = min(curReconSP);
                span   = maxVal-minVal;
                phiReconSPNorm(j,:) = (curReconSP-minVal)/span;
                
                sp(j).sigFrame = gather(A0);
                sp(j).sigPos   = gather(A5);
                sp(j).fft      = gather(A7);
                sp(j).recon    = gather(phiReconSPNorm(j,:));
                
                tAOHad = tic;
                A01 = reshape(IdHad, samplesPerSqnc, sqncPerFrame*framesPerSig);
                A02 = sMatInv * A01;
                A03 = reshape(A02, samplesPerSqnc, sqncPerFrame, framesPerSig);
                A1 = permute(A03, [2,1,3]);
                A2 = reshape(A1, sqncPerFrame, samplesPerPulse, numOfPos, framesPerSig);
                A3 = permute(A2, [2,1,3,4]);
                A4 = reshape(A3, samplesPerPos, numOfPos, framesPerSig);                                   % AC-coupling for the entire signal
                A5 = A4 - repmat(mean(A4, 1), samplesPerPos, 1, 1);       % AC-coupling for each channel - not sure needed when noise is present
                A6 = abs(fftshift(fft(A5,[],1),1)).^2;
                A7 = mean(A6, 3);
                curReconHad = gather(sqrt(A7(fUsIdxPos , :)));
                fprintf("SpeckleSim: ReconHad %.2f[s]\n", toc(tFilter));
                
                phiReconHad(j,:) = curReconHad;
                 
                maxVal = max(curReconHad);
                minVal = min(curReconHad);
                span   = maxVal-minVal;
                phiReconHadNorm(j,:) = (curReconHad-minVal)/span;
                
                had(j).sigFrame = gather(A03);
                had(j).sigPos   = gather(A5);
                had(j).fft      = gather(A7);
                had(j).recon    = gather(phiReconHadNorm(j,:));
                
                IdSPMat(:,:,j)  = gather(IdSP);
                IdHadMat(:,:,j) = gather(IdHad);
                fprintf("SpeckleSim: AOI Reconstruction in %.2f[s]\n\n", toc(trun))
                clear IdSP IdHad A01 A02 A03 A0 A1 A2 A3 A4 A5 A6 A7
            end
            
            %Collect Results:
            res.phi             = phi;
            res.phiReconHad     = phiReconHad;
            res.phiReconHadNorm = phiReconHadNorm;
            res.phiReconSP      = phiReconSP;
            res.phiReconSPNorm  = phiReconSPNorm;
            res.detSigSP        = IdSPMat;
            res.detSigHad       = IdHadMat;
            res.sp              = sp;
            res.had             = had;
            
            % collect vars
            res.vars.xUS             = this.vars.xUS;
            res.vars.fBar            = fBar;
            res.vars.fUS             = fUS;
            res.vars.fUsIdxPos       = fUsIdxPos;

            if figs
                this.displaySpeckleRes(res)
            end
        end

        function res = matchMeasAndSpeckleSim(this, phiMeas, xMeas, phiSpeckle, xSpeckle, phiConv, phi)
            
            numPh = size(phiMeas,1);
            
            figure();
            for i=1:this.vars.numMu
                ax(i) = subplot(2,3,i);
                hold on
                plot(this.vars.x, log(phi(i,:)))
                if ~isempty(phiSpeckle); plot(xSpeckle, log(abs(phiSpeckle(i,:)))); end
                if ~isempty(phiConv); plot(this.vars.x, log(abs(phiConv(i,:)))); end
                plot(xMeas, log(phiMeas(i,:)));
                xlabel("X[mm]")
                ylabel("Fluence [AU]")
                legend("MCX", "Speckle", "Conv", "Meas", 'Location', 'northwest')
                title(sprintf("Phantom: %d", i))
                xlim([-50,50])
                ylim([-10,0])
            end
            linkaxes(ax)
            
%             res.phiMeas    = phiMeas;
%             res.xMeas      = xMeas;
%             res.phi        = phi;
%             res.xPhi       = xPhi;
%             res.phiSim     = phiSim;
%             res.phiSimNorm = phiSimNorm;
%             res.xSim       = this.vars.xUS;
        end
        
        function displayResults(this)
            fprintf("VAOS: Displaying Results Figures\n");
            x = this.vars.x;
            %-------------------------------
            % Raw US Profile
            %-------------------------------
            
            figure();
            subplot(1,2,1)
            plot(this.usVars.pulseAx, this.us.pulse);
            title("Pulse Envelope")
            xlabel("X[mm]")
            ylabel("Normalized Pressure")
            subplot(1,2,2)
            plot(this.usVars.profileDepthVecRaw, this.us.focalProfileNorm);hold on
            plot(this.usVars.profileDepthVecRaw, this.us.depthProfileNorm)
            title("US Beam Axial Profile");
            xlabel("X[mm]")
            ylabel("Normalized Pressure")
            
            %-------------------------------
            % Mathematical Fluence Profile
            %-------------------------------
            figure();
            subplot(1,2,1)
            plot(x, this.phiMath')
            title("Mathematical Fluence Profile - with Reflection")
            xlabel("X[mm]")
            ylabel("Fluence [AU]")
            subplot(1,2,2)
            plot(x, log(this.phiMath'))
            title("Mathematical Fluence Profile - with Reflection")
            xlabel("X[mm]")
            ylabel("Fluence [Log]")
            
            %-------------------------------
            % Mathematical Fluence Profile
            %-------------------------------
            figure()
            hold on
            plot(x, this.phiMath);
            plot(x, this.profiles.focalProfileInt)
            plot(x, this.profiles.usIntProfile);
            xlabel("X[mm]")
            ylabel("AU")

            %------------------------------
            % Types of Pulses
            %------------------------------
            figure();
            subplot(2,2,1)
            plot(this.pulses.pulseIdeal)
            title("Ideal Delta")
            subplot(2,2,2)
            plot(this.pulses.pulseSinSig); hold on
            plot(this.pulses.pulseSin)
            title("Single Cycle Sine")
            subplot(2,2,3)
            plot(this.pulses.pulseMeas); hold on
            plot(this.pulses.pulseMeasSig)
            title("Measured Transducer")
           
            %-----------------------------
            % Spatial Single Pulse Matrix
            %-----------------------------
            figure();
            subplot(1,2,1)
            imagesc(this.spatMat.transMatEnvNormSP);
            axis tight equal
            colorbar
            subplot(1,2,2)
            imagesc(this.spatMat.transMatSigNormSP);
            axis tight equal
            colorbar

            %-----------------------------
            % Spatial Hadmarad Matrix
            %-----------------------------
            figure();
            subplot(1,2,1)
            imagesc(this.spatMat.transMatEnvNormHad)
            axis tight equal
            title("Envelope Multiplexing Mat");
            subplot(1,2,2)
            imagesc(this.spatMat.transMatSigNormHad)
            axis tight equal
            title("Signal Multiplexing Mat");
            
            %-----------------------------
            % Spatial Hadmarad Matrices
            %-----------------------------
            figure();
            subplot(2,2,1)
            imagesc(this.vars.sMat)
            title("sMat")
            colorbar
            axis tight equal
            subplot(2,2,2)
            imagesc(this.vars.sMatInv)
            axis tight equal
            title("sMatInv")
            colorbar
            subplot(2,2,3)
            imagesc(this.spatMat.transMatEnvNormHad)
            title("transMatNorm")
            axis tight equal
            colorbar
            subplot(2,2,4)
            imagesc(this.spatMat.sMatInvSqnc)
            title("sMatInvSqnc")
            axis tight equal
            colorbar
            
        end
        
        function displaySpeckleRes(this, res)
            numPh = size(res.phi,1);
            sp        = res.sp;
            had       = res.had;
            xUS       = res.vars.xUS;
            fBar      = res.vars.fBar;
            fUsIdxPos = res.vars.fUsIdxPos;
            fUS       = res.vars.fUS;
            x         = this.vars.x;
            
            
            idx = 122;
            for j=1:numPh
                figure();
                subplot(2,3,1);
                plot(squeeze(sp(j).sigFrame(:,1,1)), '-+');
                yyaxis right
                plot(squeeze(had(j).sigFrame(:,1,1)), '-o');
                title("Signal - SP"); xlabel("t[samples]"); legend("SP", "Had")
                subplot(2,3,2);
                plot(squeeze(sp(j).sigPos(:,idx,1)), '-+'); hold on
                plot(squeeze(had(j).sigPos(:,idx,1)), '-o');
                title("Signal"); xlabel("t[samples]"); legend("SP", "Had");
                subplot(2,3,3);
                plot(fBar, sp(j).fft(:, idx)); hold on
                plot(fUS*1e-6, sp(j).fft(fUsIdxPos, idx), '+g');
                yyaxis right
                plot(fBar, had(j).fft(:, idx));
                plot(fUS*1e-6, had(j).fft(fUsIdxPos, idx), '+g');
                xlabel("f[MHz]"); title("Fourier Transform"); legend("SP", "Had");
                subplot(2,3,4);
                plot(x, res.phi(j,:)); hold on
                plot(xUS, sp(j).recon, '-+');
                plot(xUS, had(j).recon, '-o');
                legend("Phi", "SP", "Had"); title("Sim. Recon"); xlabel("X[mm]")
                subplot(2,3,5);
                plot(x, log(res.phi(j,:))); hold on
                plot(xUS, log(sp(j).recon), '-+');
                plot(xUS, log(had(j).recon), '-o');
                legend("Phi", "SP", "Had"); title("Sim. Recon (Log)"); xlabel("X[mm]");
            end
            
            figure()
                for i=1:numPh
                    ax(i) = subplot(2,3,i);
                    hold (ax(i), 'on')
                    plot(ax(i), x, log(res.phi(i,:)));
                    plot(ax(i), xUS, log(res.phiReconSPNorm(i,:)), '-+');
                    plot(ax(i), xUS, log(res.phiReconHadNorm(i,:)), '-o');
                    legend(ax(i), "Phi", "SP Sim", "Had Sim")
                    title(ax(i), sprintf("Phantom: %d", i))
                    xlabel(ax(i), "X[mm]")
                    ylim(ax(i), [-10,0])
                end
                linkaxes(ax);
        end
        
        function recordConv(this)
            reconSize = this.vars.reconSize;
            focalProfileInt = this.profiles.focalProfileInt;
            usRefProfile = this.profiles.usRefProfile;
            
            [pulseEnv, pulseSig] = this.choosePulseShape();
            pulseLen = length(pulseEnv);
            
            padSize  = pulseLen-1;
            endIdx   = reconSize + padSize;
            startIdx = pulseLen;
            x = this.vars.x;
            startFigIdx = 50;
            endFigIdx = 1250;
            if this.displayDebug
                hFig1 = figure();
                ax = axes(); hold(ax, 'on');
                set(hFig1, 'color', 'white');   
                set(ax, 'color', 'white');
                h1 = plot(ax, x, this.phiMath(1,:) .* this.profiles.usRefProfile); % Fluence
                h2 = plot(ax, x, this.profiles.focalProfileInt); % US Profile
                h3 = plot(ax, x, zeros(1,reconSize)); 
                ylim(ax, [-0.05,1]); xlim(ax, [x(startFigIdx), x(endFigIdx)])
                legend("Fluence", "US Field", "US Pulse")
                xlabel("X[mm]")
                ylabel("Normalized Values [AU]")
                set(hFig1, 'Position', [288,571,560,191]);
            end

            transMatSP = zeros(reconSize);
            transMatSigSP = zeros(reconSize);
            
            pulsePos    = [1, zeros(1,reconSize-1)];
            for i = startFigIdx:endFigIdx
                curSqnc          = circshift(pulsePos, (i-1));
                curSqncSpatial   = curSqnc.*focalProfileInt.*usRefProfile;
                transmissionSP  = conv(curSqncSpatial, flip(pulseEnv), 'full');
                transMatSP(:,i) = transmissionSP(startIdx:endIdx);

                transmissionSigSP  = conv(curSqncSpatial, pulseSig, 'full');
                transMatSigSP(:,i) = transmissionSigSP(startIdx:endIdx);
            end
            
            for i = startFigIdx:endFigIdx
                if this.displayDebug
                    set(h3, 'YData', transMatSP(i,:));
                    drawnow();
                    frame(i) = getframe(hFig1);
                end 
            end
            this.spatMat.transMatEnvSP     = transMatSP;
            this.spatMat.transMatSigSP     = transMatSigSP;
            this.spatMat.transMatEnvNormSP = transMatSP/max(transMatSP(:));
            this.spatMat.transMatSigNormSP = transMatSigSP/max(transMatSigSP(:));
            
            skip = 5;
            for i = startFigIdx:skip:endFigIdx
                im = frame2im(frame(i)); 
                [imind,cm] = rgb2ind(im,256);
                if i == startFigIdx
                  imwrite(imind, cm, "spUSConv.gif", 'gif', 'DelayTime', 0, 'Loopcount', inf); 
                else 
                  imwrite(imind, cm, "spUSConv.gif", 'gif', 'DelayTime', 0, 'WriteMode', 'append'); 
                end 
            end
        end
        
    end
 end