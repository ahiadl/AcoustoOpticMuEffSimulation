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
        vars;
        pulse;

        pulses;
        phiMath;
        
        res;
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
            
            uVars.usDistFromInt = 30; %[mm]
            
            uVars.useCustomUSParams = false;
            
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
            
        %% Variables Management Function:
        function uploadUS(this, input)
            fprintf("VAOS: Loading US Field\n");
            if isstring(input) || ischar(input)
                resUS = load(input);
            elseif isstruct (input)
                resUS = input;
            else
                fprintf("VAOS: Unkown US data type.\n");
                return;
            end
            
            this.us.usTransAmpMat = resUS.data.int.trP2pNorm;
            this.us.depthProfile  = resUS.data.int.depthProfNorm;
            
            this.usVars.axialVec      = resUS.grid.int.axialVec;
            this.usVars.tr1Vec        = resUS.grid.int.tr1Vec;
            this.usVars.tr2Vec        = resUS.grid.int.tr2Vec;
            this.usVars.axialFocalIdx = resUS.grid.int.axialFocalIdx;
            this.usVars.tr1FocalIdx   = resUS.grid.int.tr1FocalIdx;
            this.usVars.tr2FocalIdx   = resUS.grid.int.tr2FocalIdx;
            this.usVars.dX            = resUS.grid.int.dAx;
            this.usVars.dtr1          = resUS.grid.int.dtr1;
            this.usVars.dtr2          = resUS.grid.int.dtr2;

            this.us.pulseEnv    = resUS.data.pulses.focalPulseEnvNorm;
            this.us.pulseSig    = resUS.data.pulses.focalPulseNorm;
            this.usVars.pulseAx = resUS.grid.pulses.axialVec;

            this.usVars.usFocalLen = resUS.usStats.focalLen;
            this.usVars.c          = resUS.usStats.c;
            
            this.vars.tVec = resUS.grid.tVec;
            this.vars.dt   = resUS.grid.dt;
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
            this.vars.usDistFromInt = uVars.usDistFromInt; %[mm] - US Focal point distance from Spacer-Phantom Interface
            
            this.vars.muEffVec = uVars.muEffVec;
            this.vars.numMu = length(this.vars.muEffVec);
        
            this.vars.useCustomUSParams = uVars.useCustomUSParams;

            this.setSpeckleVars(uVars.speckle);

            this.calcSimDimAndSpace();
            this.alignAndInterpUS();
            this.createAcousticInterface();
            this.createPulses();
        end
        
        function setSpeckleVars(this, uVars)
            this.vars.speckle.framesPerSig = uVars.framesPerSig;
            this.vars.speckle.sqncPerFrame = uVars.sqncPerFrame;
            this.vars.speckle.batchSize    = uVars.batchSize;  
            this.vars.speckle.n            = uVars.n;
            this.vars.speckle.lambda       = uVars.lambda;
            this.vars.speckle.gamma        = uVars.gamma;
            this.vars.speckle.numOfGrain   = uVars.numOfGrain;
            this.vars.speckle.SBR          = uVars.SBR;
        end

        function simVars = getVars(this)
           simVars = this.vars; 
        end
        
        %% Config Functions:
        function config(this)
            this.buildSPMatrix();
%             this.buildHadMatrix();
%             this.buildHadInvMat();
        end

        function calcSimDimAndSpace(this)
            fprintf("VAOS: Calculating Dimensions and Creating Space\n");
            
            %% Chossing US Params Source
            if this.vars.useCustomUSParams
                this.vars.us.c = 1400;
                this.vars.dX   = this.vars.dt*this.vars.us.c*1e3;
                this.vars.us.usDepthVecRaw = this.usVars.axialVec(1) + this.vars.dX*(0:1:length(this.usVars.axialVec) );
                this.vars.us.usFocalLen    = this.vars.us.usDepthVecRaw(this.usVars.axialFocalIdx);
            else
                this.vars.us.c = this.usVars.c;
                this.vars.dX   = this.usVars.dX;
                this.vars.us.usDepthVecRaw = this.usVars.axialVec;
                this.vars.us.usFocalLen    = this.usVars.usFocalLen;
            end
            
            %% Calculating Simulation Dimensions:
            this.vars.singleCycleLen    = (1/this.vars.fUS)*this.vars.us.c*1e3;
            this.vars.singleCycleIdxRaw = this.vars.singleCycleLen/this.vars.dX;
            this.vars.singleCycleIdx    = round(this.vars.singleCycleIdxRaw);
            
            if this.debug
                this.vars.sqncCompSpac = 8;
            else
                this.vars.sqncCompSpac = this.vars.singleCycleIdx;
            end
            
            this.vars.reconSize = this.vars.N*this.vars.sqncCompSpac;

            fprintf("VAOS: Creating Spatial Vars:\n");
            dX = this.vars.dX;
            reconSize = this.vars.reconSize;
            
            % Creating system axis:
            xRaw = (0:1:reconSize-1)*dX;
            x = xRaw - (this.vars.us.usFocalLen + this.vars.usDistFromInt); % mm
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
            
            spacerAirInt = floor(spacerLenIdx/2) + spacerIdxStart;
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
            
            %% Collect parameters:
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
            
        function alignAndInterpUS(this)
            fprintf("VAOS: Aligning US profile to Fluence\n");
            x = this.vars.x;
            
            [~, focalPointIdx] = min(abs(x+this.vars.usDistFromInt));
            usFocalPos = x(focalPointIdx);
            
            profileDepthVecArt = this.vars.us.usDepthVecRaw - this.vars.us.usFocalLen + usFocalPos;
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
            this.pulses.pulseMeas    = this.us.pulseEnv;
            this.pulses.pulseMeasSig = this.us.pulseSig;
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

            this.pulse.sig = sig;
            this.pulse.env = env;
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
                curSqnc            = circshift(pulsePos, (i-1));
                curSqncSpatial     = curSqnc.*focalProfileInt.*usIntProfile;
                transmissionEnvSP  = conv(curSqncSpatial, flip(pulseEnv), 'full');
                transMatEnvSP(:,i) = transmissionEnvSP(startIdx:endIdx);

                transmissionSigSP  = conv(curSqncSpatial, flip(pulseSig), 'full');
                transMatSigSP(:,i) = transmissionSigSP(startIdx:endIdx);

                if this.displayDebug && ((i==1) || (~mod(i,20)) || (i ==reconSize))
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
            
            sVecInt        = zeros(sqncCompSpac,N);
            sVecInt(1, :)  = this.vars.sVec;
            sVecInt        = logical(sVecInt(:)');
            
            % Calculate Convolution indices:
            padSize = length(pulseEnv)-1;
            transFullLen = 3*reconSize+padSize;
            startIdx = reconSize + 1;
            endIdx   = 2*reconSize;
            
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
                curIdxs = idxVec(curSqnc);
                transEnvMat(i,:) = sum(transMatEnvSP(curIdxs,:),1);
                transSigMat(i,:) = sum(transMatSigSP(curIdxs,:),1);

                if this.displayDebug  && ((i==1) || (~mod(i,20)) || (i ==reconSize))
                    set(h1, 'YData', curSqnc);
                    set(h3, 'CData', transEnvMat);
                    set(h5, 'YData', transEnvMat(i,:));
                    set(h6, 'CData', transSigMat);
                    set(h8, 'YData', transSigMat(i,:));
                    drawnow();
                end
            end
            
            this.spatMat.transMatEnvHad     = gather(transEnvMat);
            this.spatMat.transMatSigHad     = gather(transSigMat);
            this.spatMat.transMatEnvNormHad = gather(transEnvMat/ max(transEnvMat(:)));
            this.spatMat.transMatSigNormHad = gather(transSigMat/ max(transSigMat(:)));
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
        
        function usTranEnv = cutTransUSEnv (this, tr1EnvSize, tr2EnvSize)
            tr1Idx = this.usVars.tr1FocalIdx;
            tr2Idx = this.usVars.tr2FocalIdx;

            tr1Idxs = tr1Idx-tr1EnvSize:tr1Idx+tr1EnvSize;
            tr2Idxs = tr2Idx-tr2EnvSize:tr2Idx+tr2EnvSize;

            usTranEnv = this.us.usTransAmpMat(tr1Idxs,tr2Idxs);
            usTranEnv = usTranEnv/max(usTranEnv(:));

        end
        %% Simulations Functions:
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
            phiRef = exp(-muEffVec.*abs(refX-refXStart));
            
            phi = [phiInc, spacer, phiRef];

%             if this.debug
%                 reconSize = this.vars.reconSize;
%                 phi = (1e-2)*ones(numMu, reconSize);
%                 phi(:,1) = 1e-3;
%                 phi(:, 300:330) = 1;
%             end
            
            this.phiMath = phi;
        end

        function res = reconNaive(this, phi, figs)
            fprintf("VAOS: Reconstructing Naive Approach\n");
            numPh = size(phi,1);
            reconSize = this.vars.reconSize;

            tr1Size = size(phi, 3);
            tr2Size = size(phi, 4);

            tr1EnvSize = floor(tr1Size/2);
            tr2EnvSize = floor(tr2Size/2);

            phiEnvRecon     = zeros(numPh, reconSize);
            phiEnvReconNorm = zeros(numPh, reconSize);
            phiSigRecon     = zeros(numPh, reconSize);
            phiSigReconNorm = zeros(numPh, reconSize);

            trEnvMat = this.cutTransUSEnv(tr1EnvSize, tr2EnvSize);
            trEnvMat = repmat(permute(trEnvMat, [3,1,2]), reconSize, 1, 1);

            for i=1:numPh
                curPhi = squeeze(phi(i,:,:,:));
                curPhi = curPhi .* trEnvMat;
                curPhi = reshape(curPhi, reconSize, []);

                phiEnvRecon(i,:)     = sum(curPhi, 2);
                phiEnvReconNorm(i,:) = phiEnvRecon(i,:)./max(phiEnvRecon(i,:),[], 2);
            end

            res.phiEnvRecon     = phiEnvRecon;
            res.phiEnvReconNorm = phiEnvReconNorm;
            res.phiSigRecon     = phiSigRecon;
            res.phiSigReconNorm = phiSigReconNorm;

            if figs
                this.displayConvRecon(phi, phiEnvReconNorm, 'Naive Recon Env')
            end
        end

        function res = reconSP(this, phi, figs)
            fprintf("VAOS: Reconstructing with Single-Pulse US.\n")
            numPh = size(phi,1);
            reconSize = this.vars.reconSize;

            tr1Size = size(phi, 3);
            tr2Size = size(phi, 4);

            tr1EnvSize = floor(tr1Size/2);
            tr2EnvSize = floor(tr2Size/2);

            phiEnvRecon      = zeros(numPh, reconSize);
            phiEnvReconAlign = zeros(numPh, reconSize);
            phiEnvReconNorm  = zeros(numPh, reconSize);
            phiSigRecon      = zeros(numPh, reconSize);
            phiSigReconAlign = zeros(numPh, reconSize);
            phiSigReconNorm  = zeros(numPh, reconSize);

            trEnvMat = this.cutTransUSEnv(tr1EnvSize, tr2EnvSize);
            trEnvMat = repmat(permute(trEnvMat, [3,1,2]), reconSize, 1, 1);

            for i=1:numPh
                curPhi = squeeze(phi(i,:,:,:));
                curPhi = curPhi .* trEnvMat;
                curPhi = reshape(curPhi, reconSize, []);

                % Envelope Reconstruction
                phiEnvRecon(i,:)      = sum(this.spatMat.transMatEnvNormSP * curPhi, 2);
                phiEnvReconAlign(i,:) = this.interpAlignReplPhi(phiEnvRecon(i,:), this.vars.x, false, true, false);
                phiEnvReconNorm(i,:)  = phiEnvReconAlign(i,:) ./ max(phiEnvReconAlign(i,:),[], 2);
                
                % Signal Reconstruction (meaningless)
                phiSigRecon(i,:)     = sum(this.spatMat.transMatSigNormSP * curPhi, 2);
                phiSigReconAlign(i,:) = this.interpAlignReplPhi(phiSigRecon(i,:), this.vars.x, false, true, false);
                phiSigReconNorm(i,:) = phiSigRecon(i,:) ./ max(phiSigRecon(i,:),[], 2);
            end

            res.phiEnvRecon      = phiEnvRecon;
            res.phiEnvReconAlign = phiEnvReconAlign;
            res.phiEnvReconNorm  = phiEnvReconNorm;
            res.phiSigRecon      = phiSigRecon;
            res.phiSigReconAlign = phiSigReconAlign;
            res.phiSigReconNorm  = phiSigReconNorm;

            if figs
                this.displayConvRecon(phi, phiEnvReconNorm, 'SP Recon Envelope')
            end

            this.res.conv.SP = res;
        end

        function res = reconHad(this, phi, figs)
            fprintf("VAOS: Reconstructing with Hadamard US.\n")
            numPh = size(phi,1);
            reconSize = this.vars.reconSize;

            tr1Size = size(phi, 3);
            tr2Size = size(phi, 4);

            tr1EnvSize = floor(tr1Size/2);
            tr2EnvSize = floor(tr2Size/2);

            phiEnvRecon     = zeros(numPh, reconSize);
            phiEnvReconNorm = zeros(numPh, reconSize);
            phiSigRecon     = zeros(numPh, reconSize);
            phiSigReconNorm = zeros(numPh, reconSize);

            trEnvMat = this.cutTransUSEnv(tr1EnvSize, tr2EnvSize);
            trEnvMat = repmat(permute(trEnvMat, [3,1,2]), reconSize, 1, 1);
            
            hadEnvMat = gpuArray(this.spatMat.transMatEnvNormHad);
            hadSigMat = gpuArray(this.spatMat.transMatSigNormHad);
            hadInvMat = gpuArray(this.spatMat.sMatInvSqnc);

            for i=1:numPh
                curPhi = gpuArray(squeeze(phi(i,:,:,:)));
                curPhi = curPhi .* trEnvMat;
                curPhi = reshape(curPhi, reconSize, []);

                phiEnvRecon(i,:)     = gather(sum(hadInvMat * hadEnvMat * curPhi, 2));
                phiEnvReconNorm(i,:) = phiEnvRecon(i,:)./max(phiEnvRecon(i,:),[], 2);
                
                phiSigRecon(i,:)     = gather(sum(hadInvMat * hadSigMat * curPhi, 2));
                phiSigReconNorm(i,:) = phiSigRecon(i,:)./max(phiSigRecon(i,:),[], 2);
            end

            res.phiEnvRecon     = phiEnvRecon;
            res.phiEnvReconNorm = phiEnvReconNorm;
            res.phiSigRecon     = phiSigRecon;
            res.phiSigReconNorm = phiSigReconNorm;
            
            if figs
                this.displayConvRecon(phi, phiEnvReconNorm, 'Had Recon Env')
            end

            this.res.conv.had = res;
        end
        
        function res = reconAll(this, phi, depthVec, figs)
            res.fluence.phi = this.interpAlignReplPhi(phi, depthVec, true, true, false);

            res.naive   = this.reconNaive(res.fluence.phi, false);
            res.sp      = this.reconSP(res.fluence.phi, false);
            res.had     = this.reconHad(res.fluence.phi, false);
            res.speckle = this.speckleSim(res.fluence.phi, false);

            this.res = res;
            if this.figs
                this.displayAllRecon(res);
            end
        end

        function res = runSpeckleSim(this, uVars, phi, figs)
            this.setSpeckleSimVaes(uVars);
            res = this.speckleSim(phi, figs);
        end

        function res = speckleSim(this, phi, figs)
            phi = abs(phi);
            numPh = size(phi,1);
            tr1Size = size(phi,3);
            tr2Size = size(phi,4);
            tr1Env  = floor(tr1Size/2);
            tr2Env  = floor(tr2Size/2);

            %----------------------------------------
            % Calculate AO Reconstruction dimensions:
            %----------------------------------------
            reconSize       = this.vars.reconSize;
            
            framesPerSig    = this.vars.speckle.framesPerSig; % number of speckle-decorrelation frames equivalent to time-to-sample.
            sqncPerFrame    = this.vars.speckle.sqncPerFrame;
            pulsePerSqnc    = this.vars.N;
            samplesPerPulse = this.vars.singleCycleIdx; % theoretical
            
            samplesPerSqnc  = samplesPerPulse * pulsePerSqnc;
            samplesPerFrame = samplesPerSqnc  * sqncPerFrame;
            samplesPerPos   = samplesPerPulse * sqncPerFrame;
            
            numOfPos = pulsePerSqnc;
            
            
            batchSize    = this.vars.speckle.batchSize;
            numOfBatch   = ceil(framesPerSig/batchSize);
            framesPerSig = numOfBatch * batchSize;
            %----------------------------------------
            % Calculate Time-Space Grid:
            %----------------------------------------
            x = this.vars.x;
            dX = this.vars.dX;
            xMat = repmat(x, reconSize, 1, batchSize); % Time x Space
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
            % Calculate Transversal Grid
            %----------------------------------------
            tr1Idx = this.usVars.tr1FocalIdx;
            tr2Idx = this.usVars.tr2FocalIdx;
            
            usTransAmp = this.us.usTransAmpMat(tr1Idx-tr1Env:tr1Idx+tr1Env,tr2Idx-tr2Env:tr2Idx+tr2Env); 
            usTransAmp = usTransAmp/max(usTransAmp(:));
            %----------------------------------------
            % EM Parameters:
            %----------------------------------------
            transMatSigNormSP  = this.spatMat.transMatSigNormSP; % US Modulation
            transMatSigNormHad = this.spatMat.transMatSigNormHad;
            
            n      = this.vars.speckle.n;
            lambda = this.vars.speckle.lambda;
            k0     = 2*pi/lambda; 
            gamma  = this.vars.speckle.gamma; % modulation depth
            dnSP   = repmat(gpuArray(-gamma*transMatSigNormSP), 1, 1, batchSize); % IS dn linear? should find a reference
            dnHad  = repmat(gpuArray(-gamma*transMatSigNormHad), 1, 1, batchSize);
            
            numOfGrain = this.vars.speckle.numOfGrain;  % number of speckle grains;
            SBR        = this.vars.speckle.SBR; % Signal-to-background ratio -> does not affect the SNR (obviously)
            
            cleanSigSP  = zeros(samplesPerSqnc, framesPerSig, 'gpuArray');
            cleanSigHad = zeros(samplesPerSqnc, framesPerSig, 'gpuArray');
            IdSPMat    = zeros(samplesPerFrame, framesPerSig, numPh);
            IdHadMat   = zeros(samplesPerFrame, framesPerSig, numPh);
            
            %-----------------------------------
            % Create Electrical Circuit LPF:
            %-----------------------------------
            [~,  lpf]  = lowpass(zeros(1,reconSize), 2*fUS, fs, 'StopbandAttenuation', 40);
            lpfT       = lpf.Coefficients;
            lpfLen     = length(lpfT);
            ipfPivot   = floor(lpfLen/2)+1;
            padSize    = ipfPivot-1;
            paddedSize = samplesPerFrame + 2*padSize;
            
            filterMat  = zeros(samplesPerFrame, paddedSize);
            
            % Create a batch:
            batchRowSize    = 1000;
            batchColSize = batchRowSize + 2*padSize;
            batchMat     = zeros(batchRowSize, batchColSize);
            
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
            sMatInvSqnc = this.spatMat.sMatInvSqnc;
            sMatInv = gpuArray(sMatInvSqnc);
            
            %----------------------------------------
            % Init data structure
            %----------------------------------------
            phiReconSP      = zeros(numPh, numOfPos);
            phiReconHad     = zeros(numPh, numOfPos);
            phiReconSPNorm  = zeros(numPh, numOfPos);
            phiReconHadNorm = zeros(numPh, numOfPos);
            
            EgAOSP  = zeros(reconSize, 1, batchSize, 'gpuArray');
            EgAOHad = zeros(reconSize, 1, batchSize, 'gpuArray');
            
            spatPhase = gpuArray(exp(1i.*n.*k0.*xMat*1e-3));
            
            T0 = zeros(numOfBatch, tr1Size, tr2Size, 5);
            T1 = zeros(numOfBatch, tr1Size, tr2Size);
            T2 = zeros(numOfBatch, tr1Size);
            T3 = zeros(numOfBatch, 1);
            T4 = zeros(numOfBatch, 6);
            %----------------------------------------
            % Simulate:
            %----------------------------------------
            for j=1:numPh
                fprintf("SpeckleSim: Phantom #%d \n", j);
                tPh = tic;
                
                curPhi = gpuArray(phi(j,:,:,:));
%                 Ephii  = zeros(reconSize, reconSize, batchSize, 'gpuArray');
                Ebkg   = SBR*sum(sqrt(curPhi(1,:)))* ones(reconSize, numOfGrain, 'gpuArray');
                %----------------------------------------------------------
                % Intereferece between modulated and unmodulated fields:
                %----------------------------------------------------------
                tFrames = tic;
                for i=1:numOfBatch
                    
                    if i==1 || ~mod(i, 25)
                       fprintf("Batch: %d\n", i); 
                    end
                    
                    tFrame = tic;
                    
                    tAlloc = tic;
                    EgAOSP(:)  = 0;
                    EgAOHad(:) = 0;
                    T4(i,1) = toc(tAlloc);
                    
                    tEnv = tic;
                    for m = 1:tr1Size
                        tLine = tic;
                        for n = 1:tr2Size
                            tPoint = tic;
                            startIdx =  (i-1)*batchSize+1;
                            endIdx   =  i*batchSize;
                            
                            tEphi = tic;
%                             Ephii(:,:,:) = repmat(curPhi(1,:,m,n), reconSize, 1, batchSize);
                            Ephii = repmat(curPhi(1,:,m,n), reconSize, 1, batchSize);
                            T0(i,m,n,1) = toc(tEphi);
                            
                            tPhase = tic;
                            phaseX = repmat (2*pi*rand(1,numOfPos, batchSize, 'gpuArray'), samplesPerPulse, 1, 1);
                            phaseX = 1i.*repmat (reshape(phaseX, 1, samplesPerSqnc, batchSize), reconSize, 1, 1);
                            T0(i,m,n,2) = toc(tPhase);

                            tField = tic;
                            EiSP  = Ephii .*  spatPhase .* exp((1i.*usTransAmp(m,n).*k0.*dX*1e-3).*dnSP  + phaseX);
                            EiHad = Ephii .*  spatPhase .* exp((1i.*usTransAmp(m,n).*k0.*dX*1e-3).*dnHad + phaseX);
                            T0(i,m,n,3) = toc(tField);
                            
                            tGrains = tic;
                            EgAOSP(:,1,:)  = EgAOSP  + sum(EiSP,2);
                            EgAOHad(:,1,:) = EgAOHad + sum(EiHad,2);
                            T0(i,m,n,4) = toc(tGrains);
                            
                            tClear = tic;
                            clear EiHad EiSP phaseX;
                            T0(i,m,n,5) = toc(tClear);
                            
                            T1(i,m,n) = toc(tPoint);
                        end
                        T2(i,m) = toc(tLine);
                    end
                    T3(i) = toc(tEnv);

                    tBkg = tic;
                    Ebkgg   = Ebkg.*exp(repmat(1i.*2.*pi.*rand(1, numOfGrain, batchSize, 'gpuArray'), reconSize, 1, 1));
                    T4(i,2) = toc(tBkg);

                    tInt = tic;
                    EgSP  = repmat(EgAOSP,  1, numOfGrain, 1)  + Ebkgg;
                    EgHad = repmat(EgAOHad, 1, numOfGrain, 1)  + Ebkgg;
                    T4(i, 3) = toc(tInt); 

                    % Calculate net signal that hits the detector:
                    tDet = tic;
                    cleanSigSP(:, startIdx:endIdx)  = squeeze(sum(abs(EgSP.*conj(EgSP)), 2));
                    cleanSigHad(:, startIdx:endIdx) = squeeze(sum(abs(EgHad.*conj(EgHad)), 2));
                    T4(i, 4) = toc(tDet);
                    
                    tClearFrame = tic;
                    clear Ebkgg;
                    clear EgSP EgHad;
                    T4(i, 5) = toc(tClearFrame);
                    
                    T4(i, 6) = toc(tFrame);
                end
                T5 = toc(tFrames);
                
                %----------------------------------------------------------
                % Logistics calculations & clearing:
                %----------------------------------------------------------
                %Average timings
                avgT0 = 1e3*squeeze(mean(mean(mean(T0,1),2),3)); 
                avgT1 = 1e3*squeeze(mean(mean(mean(T1,1),2),3)); 
                avgT2 = mean(T2(:)); 
                avgT3 = mean(T3);
                avgT4 = 1e3*mean(T4,1);
                
                % Clear unneeded heavy variable:
                tClearPH = tic;
                
                
                clear Ephii Ebkg;
%                 clear Ebkgg;
                tClearPH = toc(tClearPH);
                
                %Print Timing:
                fprintf("SpeckleSim: Timings:\n"); 
                fprintf("(1) Average Field Calc. Timings: \n    Ephi: %.2f[ms] | Phase: %.2f[ms] | Field: %.2f[ms] | Grains: %.2f[ms] | Clear: %.2f[ms].\n",...
                        avgT0(1), avgT0(2), avgT0(3), avgT0(4), avgT0(5))
                fprintf("(2) Create Random AO Fields for:\n    Point: %.2f[ms] | Line: %.2f[s] | Env: %.2f[s].\n", avgT1, avgT2, avgT3)
                fprintf("(3) Average Frame Calc. Timings: \n   Alloc: %.2f[ms] | Bkg: %.2f[ms] | Int: %.2f[ms] | Det: %.2f[ms] | Clear: %.2f[ms] | Frame: %.2f[s].\n", ...
                        avgT4(1), avgT4(2), avgT4(3), avgT4(4), avgT4(5), 1e-3*avgT4(6));
                fprintf("(4) All frames complete randomization & interference: %.2f[s].\n", T5);
                fprintf("(5) Clearing Variables  %.2f[ms]\n", tClearPH*1e3);
                
                tSignal = tic;
                %----------------------------------------------------------
                % Replicate clean signal for Sequencing:
                %----------------------------------------------------------
                tSqnc = tic;
                cleanSigFrameSP  = repmat(cleanSigSP,  sqncPerFrame, 1);
                cleanSigFrameHad = repmat(cleanSigHad, sqncPerFrame, 1);
                tSqnc = 1e3*toc(tSqnc); fprintf("(6) Replicate clean signal for Sqnc %.2f[ms]\n", tSqnc)
                
                %----------------------------------------------------------
                % Intereferece between modulated and unmodulated fields:
                %----------------------------------------------------------
                tNoise = tic;
                noise         = sqrt(mean(cleanSigSP(:))) * randn(samplesPerFrame, framesPerSig);
                noisedSigSP   = noise + cleanSigFrameSP;
                noisedSigHad  = noise + cleanSigFrameHad;
                tNoise = 1e3*toc(tNoise); fprintf("(7) Add Noise %.2f[ms]\n", tNoise)
                
                
                % Signal after AC coupling:
                tAC = tic;
                acCoupledSP  = noisedSigSP  - mean(noisedSigSP(:));
                acCoupledHad = noisedSigHad - mean(noisedSigSP(:));
                tAC = 1e3*toc(tAC); fprintf("(8) AC Coupling %.2f[ms]\n", tAC);
                
                % Signal after LPF:
                tPad = tic;
                paddedSigSP  = [acCoupledSP(1)   * ones(padSize,framesPerSig);...
                                acCoupledSP; ...
                                acCoupledSP(end) * ones(padSize,framesPerSig)];
                paddedSigHad = [acCoupledHad(1)  * ones(padSize,framesPerSig);...
                                acCoupledHad; ...
                                acCoupledHad(end)* ones(padSize,framesPerSig)];
                fprintf("(9) Padding %.2f[ms]\n", 1e3*toc(tPad));
                
                tFilter = tic;
                IdSP  = filterMat * gather(paddedSigSP);
                IdHad = filterMat * gather(paddedSigHad);
                tFilter = 1e3*toc(tFilter); fprintf("(10) Filtering %.2f[ms]\n", tFilter)
                
                tClearSig = tic;
                clear cleanSigFrameSP cleanSigFrameHad 
                clear noise noisedSigSP noisedSigHad acCoupledSP 
                clear acCoupledHad paddedSigSP paddedSigHad
                tClearSig = 1e3*toc(tClearSig); fprintf("(11) Clear signals %.2f[ms]\n", tClearSig);
                
                tSignal = toc(tSignal); fprintf("(12) Created the measured signal in %.2f[s]\n", tSignal);
                
                tRecon = tic;
                %-----------------------------------
                % AOI Reconstruction Algorithm:
                %-----------------------------------
                % SP Reconstrcution:
                tAOSP = tic;
                A0 = reshape(IdSP, samplesPerSqnc, sqncPerFrame, framesPerSig);
                A1 = permute(A0, [2,1,3]);
                A2 = reshape(A1, sqncPerFrame, samplesPerPulse, numOfPos, framesPerSig);
                A3 = permute(A2, [2,1,3,4]);
                A4 = reshape(A3, samplesPerPos, numOfPos, framesPerSig);                                   % AC-coupling for the entire signal
                A5 = A4 - repmat(mean(A4, 1), samplesPerPos, 1, 1);       % AC-Coupling for each channel - not sure needed when noise is present
                A6 = abs(fftshift(fft(A5,[],1),1)).^2;
                A7 = mean(A6, 3);
                curReconSP = gather(sqrt(A7(fUsIdxPos , :)));
                tAOSP = 1e3*toc(tAOSP); fprintf("(13) Recon Single-Pulse: %.2f[ms]\n", tAOSP);
                
                phiReconSP(j,:)  = curReconSP;
                
                maxVal = max(curReconSP);
                minVal = min(curReconSP);
                span   = maxVal-minVal;
                phiReconSPNorm(j,:) = (curReconSP-minVal)/span;
                
                sp(j).Id       = gather(IdSP);
                sp(j).sigFrame = gather(A0);
                sp(j).sigPos   = gather(A5);
                sp(j).fft      = gather(A7);
                sp(j).recon    = gather(phiReconSPNorm(j,:));
                
                % Had Reconstrcution:
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
                tAOHad = 1e3*toc(tAOHad); fprintf("(14) Recon Hadamard: %.2f[ms]\n", tAOHad);
                
                phiReconHad(j,:) = curReconHad;
                 
                maxVal = max(curReconHad);
                minVal = min(curReconHad);
                span   = maxVal-minVal;
                phiReconHadNorm(j,:) = (curReconHad-minVal)/span;
                
                had(j).Id       = gather(IdHad);
                had(j).sigFrame = gather(A03);
                had(j).sigPos   = gather(A5);
                had(j).fft      = gather(A7);
                had(j).recon    = gather(phiReconHadNorm(j,:));
                
%                 IdSPMat(:,:,j)  = gather(IdSP);
%                 IdHadMat(:,:,j) = gather(IdHad);
                
                tRecon = toc(tRecon); fprintf("(15) Total AOI Reconstruction in %.2f[s]\n", tRecon)
                
                %-----------------------------------
                % Clear reconstruction Variable:
                %-----------------------------------
                tClearRecon = tic;
                clear IdSP IdHad A01 A02 A03 A0 A1 A2 A3 A4 A5 A6 A7
                tClearRecon = 1e3*toc(tClearRecon); fprintf("(16) Clear reconstrcution variables: %.2f[ms]\n", tClearRecon);
                
                %-----------------------------------
                tPh = toc(tPh); fprintf("(17) Complete Synthesis and Reconstruction: %.2f[s]\n\n", tPh);
            end
            
            %Collect Results:
            res.recon.phiReconHad     = phiReconHad;
            res.recon.phiReconHadNorm = phiReconHadNorm;
            res.recon.phiReconSP      = phiReconSP;
            res.recon.phiReconSPNorm  = phiReconSPNorm;
            res.AOData.sp             = sp;
            res.AOData.had            = had;

            % collect vars
            res.vars.xUS             = this.vars.xUS;
            res.vars.fBar            = fBar;
            res.vars.fUS             = fUS;
            res.vars.fUsIdxPos       = fUsIdxPos;

            this.res.resSpeckle = res;

            if figs
                this.displaySpeckleRecon(res)
            end
        end
        
        %% Display Functions:
        
        function displayResults(this)

            fprintf("VAOS: Displaying Results Figures\n");
            x = this.vars.x;
            %-------------------------------
            % Raw US Profile
            %-------------------------------
            figure();
            subplot(2,2,1)
            plot(this.usVars.pulseAx, this.pulse.env);
            title("Pulse Envelope")
            xlabel("X[mm]")
            ylabel("Normalized Pressure")
            subplot(2,2,3)
            plot(this.usVars.pulseAx, this.pulse.sig);
            title("Pulse Temporal")
            xlabel("X[mm]")
            ylabel("Normalized Pressure")
            subplot(1,2,2)
            plot(this.vars.us.usDepthVecRaw, this.us.depthProfile);
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
        
        function displaySpeckleRecon(this, res)
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
                plot(squeeze(sp(j).sigFrame(:,1,1)), '-+'); hold on
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
                plot(x, res.phi(j,:,1,1)); hold on
                plot(xUS, sp(j).recon, '-+');
                plot(xUS, had(j).recon, '-o');
                legend("Phi", "SP", "Had"); title("Sim. Recon"); xlabel("X[mm]")
                subplot(2,3,5);
                plot(x, log(res.phi(j,:,1,1))); hold on
                plot(xUS, log(sp(j).recon), '-+');
                plot(xUS, log(had(j).recon), '-o');
                legend("Phi", "SP", "Had"); title("Sim. Recon (Log)"); xlabel("X[mm]");
            end
            
            figure()
            for i=1:numPh
                ax(i) = subplot(2,3,i);
                hold (ax(i), 'on')
                plot(ax(i), x, log(res.phi(i,:,1,1)));
                plot(ax(i), xUS, log(res.phiReconSPNorm(i,:)), '-+');
                plot(ax(i), xUS, log(res.phiReconHadNorm(i,:)), '-o');
                legend(ax(i), "Phi", "SP Sim", "Had Sim")
                title(ax(i), sprintf("Phantom: %d", i))
                xlabel(ax(i), "X[mm]")
                ylim(ax(i), [-10,0])
            end
            linkaxes(ax);
        end
        
        function displayConvRecon(this, phi, phiRecon, name)

            numPh = size(phi,1);
            tr1Idx = 1; if size(phi,3) > 1; tr1Idx = floor(size(phi,3)/2); end
            tr2Idx = 1; if size(phi,4) > 1; tr2Idx = floor(size(phi,4)/2); end
            
            cols = ceil(sqrt(numPh));
            rows = floor(sqrt(numPh));

            hFig1 = figure();
            set(hFig1, 'NumberTitle', 'off', 'Name', name);
            for i=1:numPh
                ax(i) = subplot(rows,cols,i); %#ok<AGROW>
                hold on
                plot(this.vars.x, log(phi(i,:, tr1Idx, tr2Idx)))
                plot(this.vars.x, log(abs(phiRecon(i,:))) )
                xlabel("X[mm]")
                ylabel("Fluence [AU]")
                title(sprintf("Phantom: %d", i))
                xlim([-60,30]);
                ylim([-10,0]);
                legend("MCX", "Conv")
            end
            linkaxes(ax);
        end
        
        function displayAllRecon(this, res)
            if nargin < 2; res = this.res; end

            phi        = res.fluence.phi;
            naiveRes   = res.naive.phiEnvReconNorm;
            spRes      = res.sp.phiEnvReconNorm;
            
            simHad     = isfield(res, 'had'); 
            if simHad; hadRes = res.had.phiEnvReconNorm; end
            
            simSpec = isfield(res, 'speckle'); 
            if simSpec; speckleRes = res.speckle.recon.phiReconSPNorm; end

            numPh  = size(spRes, 1);
            tr1Idx = 1; if size(phi,3) > 1; tr1Idx = floor(size(phi,3)/2); end
            tr2Idx = 1; if size(phi,4) > 1; tr2Idx = floor(size(phi,4)/2); end
            
            cols = ceil(sqrt(numPh));
            rows = floor(sqrt(numPh));

            legStr = ["MCX", "Naive", "Conv SP"];
            if simHad; legStr  = [legStr, "Conv Had"]; end
            if simSpec; legStr = [legStr, "Speckle SP"]; end

            hFig1 = figure();
            set(hFig1, 'NumberTitle', 'off', 'Name', 'All VAOS Results');
            for i=1:numPh
                ax(i) = subplot(rows,cols,i); %#ok<AGROW>
                hold on
                plot(ax(i), this.vars.x, log(phi(i,:, tr1Idx, tr2Idx)) )
                plot(ax(i), this.vars.x, log(naiveRes(i,:)) )
                plot(ax(i), this.vars.x, log(abs(spRes(i,:))) );
                if simHad;  plot(ax(i), this.vars.x, log(abs(hadRes(i,:))) ); end
                if simSpec; plot(ax(i), res.speckle.vars.xUS, log(abs(speckleRes(i,:))) ); end
                xlabel("X[mm]")
                ylabel("Fluence [AU]")
                title(sprintf("Phantom: %d", i))
                legend(legStr, 'Location', 'northwest');
                xlim([-60,30]);
                ylim([-10,0]);
            end
            linkaxes(ax);
        end

        %% Misc
        function phiAligned = interpAlignReplPhi(this, phi, depthVec, rep, align, figs)
            dX = this.vars.dX;
            x = this.vars.x;
            numPhi  = size(phi, 1);
            tr1Size = size(phi, 3);
            tr2Size = size(phi, 4);
            
            spacerLen = this.vars.space.spacerLenIdx;

            %-------------------------------------
            % 1. Interpolate for x vec resolution:
            %-------------------------------------
            depthVecSimInt = depthVec(1) : dX : depthVec(end);
            phiInt1 = zeros(numPhi, length(depthVecSimInt));
            phiint1Len = length(depthVecSimInt);

            for i=1:numPhi
                for j=1:tr1Size
                    for k = 1:tr2Size
                        phiInt1(i,:,j,k) = interp1(depthVec, phi(i,:,j,k), depthVecSimInt, 'pchip');
                    end
                end
            end
            
            %-------------------------------------
            % 2. Replicate Phi to create spacer effect:
            %-------------------------------------
            if rep
                % Replicate and add spacer between replications
                phiRep        = [flip(phiInt1,2), zeros(numPhi, spacerLen, tr1Size, tr2Size), phiInt1];
                depthVec2Side = (0:1:(size(phiRep,2)-1))*dX;
                xSim          = depthVec2Side - depthVec2Side(phiint1Len);
            else
                % Don't replicate and take it as it is
                phiRep = phiInt1;
                xSim = depthVecSimInt;
            end
            
            %-------------------------------------
            % 3. Set edges of phi to zero
            %-------------------------------------
            if x(1) < xSim(1)
                phiRep = [zeros(numPhi,1, tr1Size, tr2Size), phiRep]; % this is to set the **first** point of fluence to 0;
                xSim = [x(1), xSim];
            end
            
            if x(end) > xSim(end)
                phiRep    = [phiRep, zeros(numPhi,1,tr1Size, tr2Size)];  % this is to set the **last** point of fluence to 0;
                xSim = [xSim, x(end)];
            end
            
            %-------------------------------------
            % 4. Interpolate the replicated phi with hard edges:
            %-------------------------------------
            phiInt2 = zeros(numPhi, length(x), tr1Size, tr2Size);
            for i=1:numPhi
                for j=1:tr1Size
                    for k = 1:tr2Size
                        curPhi = interp1(xSim, phiRep(i,:,j,k), x, 'pchip');
                        phiInt2(i,:,j,k) = curPhi/max(curPhi);
                    end
                end 
            end

            %-------------------------------------
            % 5. Align all peaks together:
            %-------------------------------------
            if align
                origin = find(x == 0);
                [~, I] = max(phiInt2, [], 2);
                shiftFact = origin- permute(I, [1,3,4,2]);
                phiAligned = zeros(size(phiInt2));
                for i=1:numPhi
                    for j=1:tr1Size
                        for k = 1:tr2Size
                            phiAligned(i,:,j,k) = circshift(phiInt2(i,:,j,k),  shiftFact(i,j,k));
                        end
                    end 
                end
            else
                phiAligned = phiInt2;
            end
            
            phiAligned = abs(phiAligned);
            %-------------------------------------
            % 6. Display:
            %-------------------------------------
            if figs
                figure()
                subplot(1,2,1)
                plot(x, phiAligned(:,:, floor(tr1Size/2)+1, floor(tr2Size/2)+1))
                subplot(1,2,2)
                plot(x, log(normMatf(phiAligned(:,:, floor(tr1Size/2)+1, floor(tr2Size/2)+1),2)))
            end
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
        
        function resLoaded = loadVAOSRes(this, path)
            resLoaded = load(path);
            
            this.vars     = resLoaded.vars;
            this.res      = resLoaded.res; 
            this.profiles = resLoaded.profiles;
            this.spatMat  = resLoaded.sparMat;
            this.pulse    = resLoaded.pulse;
            this.us       = resLoaded.us;
            this.usVars   = resLoaded.usVars;
        end
        
        function saveVAOSRes(this, path)
            resSave.vars     = this.vars;
            resSave.res      = this.res;
            resSave.profiles = this.profiles;
            resSave.spatMat  = this.spatMat;
            resSave.pulse    = this.pulse;
            resSave.usVars   = this.usVars;
            resSave.us       = this.us;

            timeStamp   = datetime('now','TimeZone','local','Format','d-MMM-y HH-mm-s');
            filename    = sprintf("%s/%s-VirtualAOSim-Result.mat", path, timeStamp);
            save(filename, '-Struct', 'resSave', '-v7.3');
        end
        
        function data = getData(this)
            data = this.res;
        end

    end
end

%% Drafts
%         function matchMeasAndSim(this, phiMeas, xMeas, phiConv, phi)
%             x = this.vars.x;
%             measAlign = this.interpAlignReplPhi(phiMeas, xMeas, false, true, false);
%             
%             numPh = size(phiMeas,1);
%             cIdx = 3;
%             legStr = ["MCX", "Conv", "Meas"];
%             
%             if isempty(phiSpeckle); legStr(2) = []; cIdx =2; end
%             if isempty(phiConv);   legStr(cIdx) = []; end
%             
%             figure();
%             for i=1:this.vars.numMu
%                 ax(i) = subplot(2,3,i);
%                 hold on
%                 plot(x, log(phi(i,:)))
%                 if ~isempty(phiSpeckle); plot(x, log(abs(speckSPAlign(i,:)))); end
%                 if ~isempty(phiConv); plot(x, log(abs(phiConv(i,:)))); end
%                 plot(xMeas, log(phiMeas(i,:)));
%                 xlabel("X[mm]")
%                 ylabel("Fluence [AU]")
%                 legend(legStr, 'Location', 'northwest')
%                 title(sprintf("Phantom: %d", i))
%                 xlim([-50,50])
%                 ylim([-10,0])
%             end
%             linkaxes(ax)
% 
%         end
        

%         function res = matchMeasAndSpeckleSim(this, phiMeas, xMeas, phiSpeckle, xSpeckle, phiConv, phi)           
%             x = this.vars.x;
%             
%             measAlign = this.interpAlignReplPhi(phiMeas, xMeas, false, true, false);
%             if ~isempty(phiSpeckle)
%                 speckSPAlign = this.interpAlignReplPhi(phiSpeckle, xSpeckle, false, true, false);
%             end
%             
%             numPh = size(phiMeas,1);
%             cIdx = 3;
%             legStr = ["MCX", "Speckle", "Conv", "Meas"];
%             
%             if isempty(phiSpeckle); legStr(2) = []; cIdx =2; end
%             if isempty(phiConv);   legStr(cIdx) = []; end
%             
%             figure();
%             for i=1:this.vars.numMu
%                 ax(i) = subplot(2,3,i);
%                 hold on
%                 plot(x, log(phi(i,:)))
%                 if ~isempty(phiSpeckle); plot(x, log(abs(speckSPAlign(i,:)))); end
%                 if ~isempty(phiConv); plot(x, log(abs(phiConv(i,:)))); end
%                 plot(xMeas, log(phiMeas(i,:)));
%                 xlabel("X[mm]")
%                 ylabel("Fluence [AU]")
%                 legend(legStr, 'Location', 'northwest')
%                 title(sprintf("Phantom: %d", i))
%                 xlim([-50,50])
%                 ylim([-10,0])
%             end
%             linkaxes(ax)
%             
%             res.measAlign = measAlign;
%         end