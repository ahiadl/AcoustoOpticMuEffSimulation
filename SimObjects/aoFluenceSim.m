classdef aoFluenceSim <handle
    % Formula for muEff from:
    % [1] "Biomedical Optics Principle and Imaging", L.V.Wang, 2008, pages 98-99 (108-109).

    properties
        usa
        ia
        fs
        vaos
        ml

        data
        vars
        grid
        paths
        oprtn
        
        rawDB

    end
    
    methods (Static)
        function uVars = createUserVars()
            % Save:
            uVars.save.savePath = '.';
            uVars.save.saveFlag = false;
            uVars.save.simName  = [];
            
            % Global:
            uVars.muaVec = []; % mm^-1
            uVars.musVec = []; % mm^-1
            uVars.g      = []; 
            uVars.ref    = []; 

            % Fluence:
            uVars.fluence.loadPhi       = false;
            uVars.fluence.phiPath       = [];
            uVars.fluence.srcPath       = [];
            uVars.fluence.detPath       = [];
            uVars.fluence.LightGeometry = 'Measured';
            uVars.fluence.srcParam      = [];
            uVars.fluence.nphoton       = 1e6; %1e8
            uVars.fluence.meshSize      = []; % lenX, lenY, lenZ
            uVars.fluence.meshRes       = 0.25; % [mm]
            
            % VAOS:
            uVars.vaos.loadVaos       = false;
            uVars.vaos.vaosPath       = [];
            uVars.vaos.usPath         = [];
            uVars.vaos.usFocalDist    = 0; %distance of the US focus from illumination plane
            uVars.vaos.N              = 251;
            uVars.vaos.fUS            = 1.25e6;
            uVars.vaos.pulseType      = 'measured';
            uVars.vaos.spacerLen      = 6.4; %mm
            uVars.vaos.spacerMaterial = 'PDMS';
            uVars.vaos.mathFluence    = false;
            uVars.vaos.envSize        = 0;

            % Measurement Loader:
            uVars.ml.measPath     = [];
            uVars.ml.loadMeas     = false;
            uVars.ml.measNameTemp = "Phantom-%d";
        
        end
    end
    
    methods
        function this = aoFluenceSim(type)
            this.fs   = fluenceSim(type);
            this.ia   = illuminationAnalysis();
            this.vaos = virtualAOSim();
            this.ml   = measLoader();
            
            this.paths.curSrcPath  = '';
            this.paths.curDetPath  = '';
            this.paths.curUSPath   = '';
            this.paths.curPhiPath  = '';
            this.paths.curVaosPath = '';
            this.paths.curMeasPath = '';

            

            this.vars.table.varNames = ["Sample No. [#]"; "Sim. MuEff [1/mm]"; "Calc. Gradient"; "Calc. MuEff[1/mm]"];
        end
        
        %% Manage Variables
        function setVars(this, uVars)
            this.setSaveVars(uVars)
            
            this.setOpticalProperties(uVars);

            this.oprtn.loadPhi  = uVars.vaos.mathFluence;
            this.paths.src      = uVars.fluence.srcPath;
            this.paths.det      = uVars.fluence.detPath;
            this.paths.phiPath  = uVars.fluence.phiPath;
            
            this.vars.muEffLevels = uVars.muEffLevels;
            
            this.vars.trEnvSize = uVars.trEnvSize;

            this.loadPhi();

            this.setVaosVars(uVars);
            this.setMLVars(uVars);
        end
        
        function setOpticalProperties(this, uVars)
            this.vars.op.muaVec = uVars.muaVec;
            this.vars.op.musVec = uVars.musVec;
            this.vars.op.g      = uVars.g;
            this.vars.op.ref    = uVars.ref;
            
            this.vars.op.musPVec = this.vars.op.musVec*(1-this.vars.op.g);
            this.vars.op.muEffGT =  sqrt(3*this.vars.op.muaVec .* (this.vars.op.muaVec + this.vars.op.musPVec));

            this.vars.numOfPhantoms = size(this.vars.op.muaVec, 2);
        end

        function setFluenceVars(this, uVars)
            fprintf("AOSim: Setting vars to fluence Sim");
            this.vars.uVars.fluence = uVars.fluence;
            this.oprtn.loadPhi  = uVars.fluence.loadPhi;
 
            this.paths.src      = uVars.fluence.srcPath;
            this.paths.det      = uVars.fluence.detPath;
            this.paths.phiPath  = uVars.fluence.phiPath;
            
            this.vars.trEnvSize = uVars.trEnvSize;

            this.setOpticalProperties(uVars);
            % Fluence Sim
            if this.oprtn.loadPhi
                fprintf("-> Loading Fluence")
                this.loadPhi();
                fprintf("-> Done!")
            else
                % Preliminary Configuration of fluence sim object to obtain
                % parameters:
                phiSimVars = fluenceSim.createUserVars();
    
                phiSimVars.mua   = this.vars.op.muaVec(1);
                phiSimVars.mus   = this.vars.op.musVec(1);
                phiSimVars.g     = this.vars.op.g;
                phiSimVars.ref   = this.vars.op.ref;
               
                phiSimVars.nphoton   = uVars.fluence.nphoton;
                
                phiSimVars.meshSize  = uVars.fluence.meshSize; % lenX, lenY, lenZ
                phiSimVars.meshRes   = uVars.fluence.meshRes; % [mm]
                
                %Static Configurations:
                phiSimVars.plane       = 'YZ';
                phiSimVars.srcDir      = [1,0,0];
                phiSimVars.loadPattern = false; 
                phiSimVars.intFactor   = 1;

                % Retreive Simulation Variables:
                this.fs.setVars(phiSimVars)
                this.fs.config();
                this.vars.fluence.uVars = phiSimVars;
                this.vars.fluence.LightGeometry = uVars.fluence.LightGeometry;
                this.vars.fluence.srcParam      = uVars.fluence.srcParam;

                this.loadOptics();
            end
            fprintf("\n");
        end

        function setVaosVars(this, uVars)
            fprintf("AOFSim: Setting vars to VAOS\n");
            this.vars.uVars.vaos = uVars.vaos;

            this.oprtn.loadVaos   = uVars.vaos.loadVaos;
            this.oprtn.simSpeckle = uVars.vaos.simSpeckle;
            this.paths.usPath     = uVars.vaos.usPath;
            this.paths.vaosPath   = uVars.vaos.vaosPath;

            if uVars.vaos.loadVaos
                fprintf("-> Loading VAOS Data")
                this.loadVAOS();
                fprintf("-> Done!")
            else

                if uVars.vaos.mathFluence
                    this.setOpticalProperties(uVars)
                    this.generateMathematicalFluence();
                end

                % Load US Databases
                this.loadUS();

                % Collect VAOS Variables:
                uVarsVAOS = virtualAOSim.createUserVars();
    
                uVarsVAOS.N              = uVars.vaos.N;
                uVarsVAOS.fUS            = uVars.vaos.fUS;
                uVarsVAOS.pulseType      = uVars.vaos.pulseType; %delta, measured, 
                uVarsVAOS.spacerLen      = uVars.vaos.spacerLen; %mm
                uVarsVAOS.spacerMaterial = uVars.vaos.spacerMaterial;
                uVarsVAOS.usDistFromInt  = uVars.vaos.usFocalDist; %[mm] 
                
                uVarsVAOS.speckle.framesPerSig = uVars.vaos.speckle.framesPerSig;
                uVarsVAOS.speckle.sqncPerFrame = uVars.vaos.speckle.sqncPerFrame;
                uVarsVAOS.speckle.batchSize    = uVars.vaos.speckle.batchSize;
                uVarsVAOS.speckle.n            = this.vars.op.ref;
                uVarsVAOS.speckle.lambda       = uVars.vaos.speckle.lambda;
                uVarsVAOS.speckle.gamma        = uVars.vaos.speckle.gamma;
                uVarsVAOS.speckle.numOfGrain   = uVars.vaos.speckle.numOfGrain;
                uVarsVAOS.speckle.SBR          = uVars.vaos.speckle.SBR;
                
                % Static debug-disabled configurations:
                uVarsVAOS.debug             = false;
                uVarsVAOS.displayDebug      = false;
                uVarsVAOS.debugTime         = false;
                uVarsVAOS.useCustomUSParams = false;
                uVarsVAOS.muEffVec          = [];

                this.vaos.setVars(uVarsVAOS);
                this.vars.vaos             = this.vaos.getVars();
                this.vars.vaos.uVars       = uVarsVAOS;
                this.vars.vaos.mathFluence = uVars.vaos.mathFluence;
                this.vars.vaos.trEnvSize   = this.vars.trEnvSize;
            end
            fprintf("\n");
        end

        function setMLVars(this, uVars)
            fprintf("AOFSim: Setting vars to Measurement Loader");
            this.vars.uVars.ml = uVars.ml;
            this.oprtn.loadMeas = uVars.ml.loadMeas;
            this.paths.measPath =  uVars.ml.measPath;

            if this.oprtn.loadMeas
                uVarsML = this.ml.createUserVars();
                
                uVarsML.projectPath  = this.paths.measPath;
                uVarsML.measNameTemp = uVars.ml.measNameTemp;
                uVarsML.numMeas      = this.vars.numOfPhantoms;
                uVarsML.loadNew      = uVars.ml.loadNew;

                uVarsML.idxLow        = uVars.ml.idxLow;
                uVarsML.idxHigh       = uVars.ml.idxHigh;
                uVarsML.phiHighToLow  = uVars.ml.phiHighToLow;
                uVarsML.noiseIdxs     = uVars.ml.noiseIdxs;

                this.ml.setVars(uVarsML);
                fprintf("-> Loading ML Data")
                this.loadMeas();
                fprintf("-> Done!")

                this.vars.ml.uVars = uVarsML;
                this.vars.ml.c     = uVars.ml.c;
                this.vars.ml.alignToSig = uVars.ml.alignToSig;
            end
            fprintf("\n");
        end
        
        function setSaveVars(this, uVars)
            % Save:
            this.vars.save.path    = uVars.savePath;
            this.vars.save.flag    = uVars.saveFlag;
            this.vars.save.simName = uVars.simName;
        end

        function vars = getVars(this)
           vars = this.vars;
        end

         %% Load Functions
        function loadOptics(this)
            % Illumination and Detection
            if ~strcmp(this.paths.curSrcPath, this.paths.src)
                fprintf("AOSim: Loading Source Pattern\n");
                [srcPattern, srcVars] = this.ia.analyse(this.paths.src);
                this.rawDB.optical.srcPattern = srcPattern;
                this.rawDB.optical.srcVars = srcVars;
                this.paths.curSrcPath = this.paths.src;
            else
                srcPattern = this.rawDB.optical.srcPattern;
                srcVars    = this.rawDB.optical.srcVars;
            end
            this.data.optics.srcPattern = srcPattern.mask;
            this.data.optics.srcPatternVec1 = srcVars.scanVec-min(srcVars.scanVec);
            this.data.optics.srcPatternVec2 = srcVars.disc1Vec-min(srcVars.disc1Vec);

            
            if ~strcmp(this.paths.curDetPath, this.paths.det)
                fprintf("AOSim: Loading Detector Pattern\n");
                [detPattern, detVars]=this.ia.analyse(this.paths.det);
                this.rawDB.optical.detPattern = detPattern;
                this.rawDB.optical.detVars = detVars;
                this.paths.curDetPath = this.paths.det;
            else
                detPattern = this.rawDB.optical.detPattern;
                detVars    = this.rawDB.optical.detVars;
            end
            this.data.optics.detPattern = detPattern.mask;
            this.data.optics.detPatternVec1 = detVars.scanVec-min(detVars.scanVec);
            this.data.optics.detPatternVec2 = detVars.disc1Vec-min(detVars.disc1Vec);

            this.ia.flush();
        end

        function loadUS(this)
            if ~isempty(this.paths.usPath) && ~strcmp(this.paths.curUSPath, this.paths.usPath) 
                fprintf("AOSim: Loading US Data\n");
                dataUS = load(this.paths.usPath);
                this.rawDB.us = dataUS;
                this.paths.curUSPath = this.paths.usPath;
            else
                dataUS = this.rawDB.us;
            end
            
            this.data.us.data  = dataUS.data;
            this.vars.us.stats = dataUS.usStats;
            this.vars.us.grid  = dataUS.grid;
            
            this.vaos.uploadUS(dataUS);
        end
        
        function loadPhi(this)
            if ~isempty(this.paths.phiPath) && ~strcmp(this.paths.curPhiPath, this.paths.phiPath)
                fprintf("AOSim: Loading Fluence Data\n");
                lf = load(this.paths.phiPath);
                this.rawDB.fluence = lf;
                this.paths.curPhiPath = this.paths.phiPath;
            else
                lf = this.rawDB.fluence;
            end

            this.data.fluence = lf.data;

            this.vars.op.muaVec = lf.vars.op.muaVec;
            this.vars.op.musVec = lf.vars.op.musVec;
            this.vars.op.g      = lf.vars.op.g;
            this.vars.op.ref    = lf.vars.op.ref;
            
            this.setOpticalProperties(lf.vars.op);
            
            this.vars.fluence = lf.vars;
            this.vars.fluence = rmfield(this.vars.fluence, 'op');
            this.vars.fluence = rmfield(this.vars.fluence, 'numOfPhantoms');

            this.extractFluenceSim();
        end
        
        function loadVAOS(this)
            if ~isempty(this.vars.vaosPath) && ~strcmp(this.paths.curVaosPath, this.paths.vaosPath) 
                this.vaos.loadVaos(this.paths.curVaosPath)
                
                this.vars.vaos = this.vaos.getVars();
                this.data.vaos = this.vaos.getData();

                this.paths.curVaosPath = this.paths.vaosPath;
            end
        end

        function loadMeas(this)
            if  ~isempty(this.paths.measPath) && ~strcmp(this.paths.curMeasPath, this.paths.measPath)
                this.ml.loadMeas();
%                 this.paths.curMeasPath = this.paths.measPath;
            end
            this.paths.curMeasPath = this.paths.measPath;
        end

        %% Processing Functions:
        function config(this)
            this.vaos.config();
            this.extractVAOSGrid();
            if this.oprtn.loadMeas
                this.analyseMeas();
            end
        end

        function extractFluenceSim(this)
            this.grid.fluence.depthVec = this.vars.fluence.srcVars.xVec;
            this.grid.fluence.tr1Vec   = this.vars.fluence.srcVars.yVec;
            this.grid.fluence.tr2Vec   = this.vars.fluence.srcVars.zVec;

            this.grid.fluence.depthSize = length(this.grid.fluence.depthVec);
            this.grid.fluence.tr1Size   = length(this.grid.fluence.tr1Vec);
            this.grid.fluence.tr2Size   = length(this.grid.fluence.tr2Vec);

            this.grid.fluence.depthMid = ceil(this.grid.fluence.depthSize /2);
            this.grid.fluence.tr1Mid   = ceil(this.grid.fluence.tr1Size /2);
            this.grid.fluence.tr2Mid   = ceil(this.grid.fluence.tr2Size /2);

            this.grid.fluence.dDepth = this.vars.fluence.srcVars.meshRes;
            this.grid.fluence.dtr1   = this.vars.fluence.srcVars.meshRes;
            this.grid.fluence.dtr2   = this.vars.fluence.srcVars.meshRes;
            
            % Cut Environment
            % TODO: think how to match transverse resolution of fluence sim
            % with US field.
            tr1Mid = this.grid.fluence.tr1Mid;
            tr2Mid = this.grid.fluence.tr2Mid;

            envSize = this.vars.trEnvSize;
           
            this.grid.fluence.tr1EnvMid = envSize+1;
            this.grid.fluence.tr2EnvMid = envSize+1;
            
            tr1EnvMid = this.grid.fluence.tr1EnvMid;
            tr2EnvMid = this.grid.fluence.tr2EnvMid;

            tr1EnvIdxs = tr1Mid-envSize:tr1Mid+envSize;
            tr2EnvIdxs = tr2Mid-envSize:tr2Mid+envSize;
        
            tr1EnvSize = length(tr1EnvIdxs);
            tr2EnvSize = length(tr2EnvIdxs);

            this.data.fluence.phiEnv     = zeros(this.vars.numOfPhantoms,...
                                             this.grid.fluence.depthSize,...
                                             tr1EnvSize,...
                                             tr2EnvSize);
            this.data.fluence.phiEnvRoot = zeros(this.vars.numOfPhantoms,...
                                             this.grid.fluence.depthSize,...
                                             tr1EnvSize,...
                                             tr2EnvSize);
            for i=1:this.vars.numOfPhantoms
                % Phi Src 
                this.data.fluence.phiSrcEnv(i, :, :, :) =...
                    this.data.fluence.phiSrcNorm{i}(:, tr1EnvIdxs, tr2EnvIdxs);
                % Phi Det
                this.data.fluence.phiDetEnv(i, :, :, :) =...
                    this.data.fluence.phiDetNorm{i}(:, tr1EnvIdxs, tr2EnvIdxs);
                % Phi Light
                this.data.fluence.phiEnv(i, :, :, :) =...
                    this.data.fluence.phiLightNorm{i}(:, tr1EnvIdxs, tr2EnvIdxs);
                % Phi Light Root
                this.data.fluence.phiEnvRoot(i, :, :, :) =...
                    this.data.fluence.phiLightNormRoot{i}(:, tr1EnvIdxs, tr2EnvIdxs);
            end
            
            this.data.fluence.phiSrcEnvLog  = log(this.data.fluence.phiSrcEnv);
            this.data.fluence.phiDetEnvLog  = log(this.data.fluence.phiDetEnv);
            this.data.fluence.phiEnvLog     = log(this.data.fluence.phiEnv);
            this.data.fluence.phiEnvRootLog = log(this.data.fluence.phiEnvRoot);

            this.data.fluence.phiSrcMidLog  = this.data.fluence.phiSrcEnvLog(:, :, tr1EnvMid, tr2EnvMid);
            this.data.fluence.phiDetMidLog  = this.data.fluence.phiDetEnvLog(:, :, tr1EnvMid, tr2EnvMid);
            this.data.fluence.phiMidLog     = this.data.fluence.phiEnvLog(:, :, tr1EnvMid, tr2EnvMid);
            this.data.fluence.phiMidRootLog = this.data.fluence.phiEnvRootLog(:, :, tr1EnvMid, tr2EnvMid);

            this.grid.fluence.tr1EnvIdxs = tr1EnvIdxs; 
            this.grid.fluence.tr2EnvIdxs = tr2EnvIdxs; 
            this.grid.fluence.tr1EnvSize = tr1EnvSize; 
            this.grid.fluence.tr2EnvSize = tr2EnvSize;
        end
        
        function muEffTotal = extractMuEffFluence(this)
            numPh = this.vars.numOfPhantoms;
            muEffGT = this.vars.op.muEffGT;
    
            name = ["src", "det", "light", "lightRoot"];
            
            dataPhi.src        =this.data.fluence.phiSrcMidLog;
            dataPhi.det        =this.data.fluence.phiDetMidLog; 
            dataPhi.light      =this.data.fluence.phiMidLog;    
            dataPhi.lightRoot  =this.data.fluence.phiMidRootLog;
            
            depthVec = this.grid.fluence.depthVec;
            dx = depthVec(2) - depthVec(1);
            
            itr = length(name);
            for j=1:numPh
                for i=1:itr
                    curName     = name(i);
                    curPhiLog   = dataPhi.(curName)(j,:);

                    % Average Gradient:
                    gradVals     = gradient(curPhiLog);
                    muEffGrad.(curName)(j) = mean(abs(gradVals(2:end-1)))/dx;
                    
                    % Fitted Log Slope:
                    if strcmp(curName, 'light')
                        fitLogRes   = fit(depthVec(20:end)', curPhiLog(20:end)'/2, 'poly1');
                    else
                        fitLogRes   = fit(depthVec(20:end)', curPhiLog(20:end)', 'poly1');
                    end
                    muEffLogFit.(curName)(j) = abs(fitLogRes.p1);
                end
            end
            muEffTotal.grad = muEffGrad;
            muEffTotal.fit  = muEffLogFit;
            
            this.data.muEff.fluence = muEffTotal;

            colNames = ["Phantom idx", "GT", "Src", "Det", "Phi", "PhiRoot"];
            fprintf("Fit Calculation:\n")
            TFit = table((1:numPh)', muEffGT', muEffLogFit.src', muEffLogFit.det', muEffLogFit.light', muEffLogFit.lightRoot', 'VariableNames', colNames);
            disp(TFit);
            
            % Calc Ratio:
            for i=1:numPh
                ratio.gt(i)    = muEffGT(i)/muEffGT(1);
                ratio.src(i)   = muEffLogFit.src(i)/muEffLogFit.src(1);
                ratio.det(i)   = muEffLogFit.det(i)/muEffLogFit.det(1);
                ratio.light(i) = muEffLogFit.light(i)/muEffLogFit.light(1);
            end
        
            fprintf("Ratio Calculation:\n");
            colNames = ["Phantom idx", "GT", "Src", "Det", "Phi"];
            TRatio = table((1:numPh)', ratio.gt', ratio.src', ratio.det', ratio.light');
            TRatio.Properties.VariableNames = colNames;
            disp(TRatio);
        end

        function muEffComp = compareMuEff(this)
            isSpec  = this.oprtn.simSpeckle;
            isMeas  = this.oprtn.loadMeas;
            numPh   = this.vars.numOfPhantoms;
            muEffGT = this.vars.op.muEffGT;
    
            levelsLog.fluence = this.vars.muEffLevels.sim;
            levelsLog.naive   = this.vars.muEffLevels.sim;
            levelsLog.conv    = this.vars.muEffLevels.sim;
            
            name     = ["fluence", "naive", "conv"];
            colNames = ["Phantom idx", "GT", "Fluence", "Naive", "Conv"];
            
            dataLog.fluence = this.data.vaos.fluence.phiMidLog;
            dataLog.naive   = this.data.vaos.naive.phiEnvLog;
            dataLog.conv    = this.data.vaos.sp.phiEnvLog;

            dataLin.fluence = this.data.vaos.fluence.phiMid;
            dataLin.naive   = this.data.vaos.naive.phiEnvReconNorm;
            dataLin.conv    = this.data.vaos.sp.phiEnvReconNorm;
            
            depthVec.fluence = this.grid.vaos.depthVec;
            depthVec.naive   = this.grid.vaos.depthVec;
            depthVec.conv    = this.grid.vaos.depthVec; 

            if isSpec
                dataLog.speckle      = this.data.vaos.speckle.recon.phiSPLog;
                dataLin.speckle      = this.data.vaos.speckle.recon.phiSPMidNorm;
                depthVec.speckle     = this.data.vaos.speckle.vars.xUS;
                levelsLog.speckle      = this.vars.muEffLevels.sim;
                name                 = [name, "speckle"];
                colNames             = [colNames, "Speckle"];
            end

            if isMeas
                dataLog.meas  = this.data.meas.phiLog;
                dataLin.meas  = this.data.meas.phiNorm;
                depthVec.meas = this.grid.meas.depthVec;
                levelsLog.meas  = this.vars.muEffLevels.meas;
                name          = [name, "meas"];
                colNames      = [colNames, "Meas"];
            end
            
            phantomLen = this.grid.fluence.depthVec(end);
            
            itr = length(name);

            % Calc For Phi
            for j=1:numPh
                for i=1:itr
                    curName      = name(i);
                    curPhiLog    = dataLog.(curName)(j,:);
                    curPhiLin    = dataLin.(curName)(j,:);
                    fullDepthVec = depthVec.(curName);
                    curLevels    = levelsLog.(curName);
                    
                    % Calculate the noise level:
                    fullIdxs   = 1:length(fullDepthVec);
                    [~, phantomStartIdx] = min(abs(fullDepthVec + phantomLen));
                    tailIdxs    = 1:(phantomStartIdx-1);
                    tailPhi     = curPhiLin(tailIdxs);
                    tailPhi(tailPhi == Inf) = [];
                    avgTailRaw = mean(tailPhi);
                    tailPhi(tailPhi > 3*avgTailRaw) = avgTailRaw;
                    noiseAvg   = mean(tailPhi);
                    noiseStd   = std(tailPhi);
                    noiseLevel = log(noiseAvg + 3*noiseStd);
                    
                    % Preliminary Fluence cur accirding to noise level:
                    tmpPhi  = curPhiLog(fullDepthVec<0);
                    curIdxs = find(tmpPhi>noiseLevel & tmpPhi<curLevels(2));

                    % Remove artifacts in tail:
                    curIdxs = setdiff(curIdxs,tailIdxs);
                    
                    % Cut The fluence Again:
                    phiLogCut      = curPhiLog(curIdxs);
                    depthVecCutLog = fullDepthVec(curIdxs);
                    
                    % Find noise Level and remove everything beyond it:
                    noiseIdxs   = find(phiLogCut<noiseLevel);
                    maxNoiseIdx = max(noiseIdxs);
                    if~isempty(maxNoiseIdx)
                        curIdxs(1:maxNoiseIdx) = [];
                    end
                    
                    phiLogCut      = curPhiLog(curIdxs);
                    depthVecCutLog = fullDepthVec(curIdxs);

                    dx = fullDepthVec(2) - fullDepthVec(1);

                    % Average Gradient:
                    gradVals     = gradient(phiLogCut/2);
                    muEffGrad.(curName)(j) = mean(abs(gradVals(2:end-1)))/dx;
                    
                    % Fitted Log Slope:
                    [fitLogRes, gof]        = fit(depthVecCutLog', phiLogCut'/2, 'poly1');
                    rsqr.(curName)(j)       = gof.rsquare;
                    logFit.(curName)(j)     = abs(fitLogRes.p1);
                    [fitLogRes, ~]          = fit(depthVecCutLog', phiLogCut', 'poly1');
                    fitModel.(curName){j,:} = fitLogRes.p1*depthVecCutLog + fitLogRes.p2;
      
                    muEffData(j).(curName).phiLog      = curPhiLog;      %#ok<AGROW> 
                    muEffData(j).(curName).depthVec    = fullDepthVec;   %#ok<AGROW> 
                    muEffData(j).(curName).phiLogCut   = phiLogCut;      %#ok<AGROW> 
                    muEffData(j).(curName).depthVecCut = depthVecCutLog; %#ok<AGROW> 
                    muEffData(j).(curName).idxs        = curIdxs;        %#ok<AGROW>
                end
            end

%             strcmp(curName, 'meas')
            muEffComp.data       = muEffData;
            muEffComp.grad       = muEffGrad;
            muEffComp.fit        = logFit;
            muEffComp.model      = fitModel;
            muEffComp.rsqr       = rsqr;
            
%             % DEBUG:
%             figure();
%             plot(fullDepthVec, curPhiLog); hold on;
%             plot(depthVecCutLog, phiLogCut)
%             plot(depthVecCutLog, fitModel.(curName){j,:})

            this.data.muEff.vaos = muEffComp;

            fprintf("Fit Calculation:\n");
            TFit = table((1:numPh)', muEffGT', logFit.fluence', logFit.naive', logFit.conv');
            if isSpec; TFit.speckle = logFit.speckle'; end
            if isMeas; TFit.meas    = logFit.meas'; end
            TFit.Properties.VariableNames = colNames;
            disp(TFit);
            
            % Calc Ratio:
            for i=1:numPh
                ratio.gt(i)      = muEffGT(i)/muEffGT(1);
                ratio.fluence(i) = logFit.fluence(i)/logFit.fluence(1);
                ratio.naive(i)   = logFit.naive(i)/logFit.naive(1);
                ratio.conv(i)    = logFit.conv(i)/logFit.conv(1);
                if isSpec
                    ratio.speckle(i) = logFit.speckle(i)/logFit.speckle(1);
                end
                if isMeas
                    ratio.meas(i) = logFit.meas(i)/logFit.meas(1);
                end
            end
        
            fprintf("Ratio Calculation:\n");
            TRatio = table((1:numPh)', ratio.gt', ratio.fluence', ratio.naive', ratio.conv');
            if isSpec; TRatio.speckle = ratio.speckle'; end
            if isMeas; TRatio.meas    = ratio.meas'; end
            TRatio.Properties.VariableNames = colNames;
            disp(TRatio);
            
            muEffComp.ratio      = ratio;
            TRatio.TFit          = TFit;
            muEffComp.TRatio     = TRatio;

            this.data.muEff.vaos = muEffComp;
        end

        function extractVAOSGrid(this)
            this.grid.vaos.depthVec   = this.vars.vaos.x;
            this.grid.vaos.dDepth     = this.vars.vaos.dX;
            this.grid.vaos.phiSize    = this.vars.vaos.reconSize;
            this.grid.vaos.phiIdxLow  = 1;
            this.grid.vaos.phiIdxHigh = this.vars.vaos.space.incIdxEnd;
        end

        function analyseMeas(this)
            % Analyse the measurement:
            this.ml.resetCurData();
            this.ml.fixSpeedOfSound(this.vars.ml.c);
            this.ml.interp(this.grid.vaos.dDepth)
            
            if this.vars.ml.alignToSig
                this.ml.alignToSignalMax();
                this.data.meas = this.ml.data.alignSig;
                this.grid.meas.depthVec = this.ml.grid.depthVecAlignSig;
                this.grid.meas.dDepth   = this.ml.grid.dDepthAlignSig;
            else
                this.ml.alignDepthVec();
                this.data.meas = this.ml.data.align;
                this.grid.meas.depthVec = this.ml.grid.depthVecAligned;
                this.grid.meas.dDepth   = this.ml.grid.dDepthAligned;
            end
%             this.ml.displayResults();          
        end

        function analyseResults(this)
            this.data.vaos.fluence.phiLog    = log(abs(this.data.vaos.fluence.phi));
            this.data.vaos.fluence.phiMid    = this.data.vaos.fluence.phi(:,:,this.grid.fluence.tr1EnvMid, this.grid.fluence.tr2EnvMid);
            this.data.vaos.fluence.phiMidLog = this.data.vaos.fluence.phiLog(:,:,this.grid.fluence.tr1EnvMid, this.grid.fluence.tr2EnvMid);

%             this.data.vaos.fluence.phiRootLog    = log(abs(this.data.vaos.fluence.phiRoot));
%             this.data.vaos.fluence.phiRootMid    = this.data.vaos.fluence.phiRoot(:,:,this.grid.fluence.tr1EnvMid, this.grid.fluence.tr2EnvMid);
%             this.data.vaos.fluence.phiRootMidLog = this.data.vaos.fluence.phiRootLog(:,:,this.grid.fluence.tr1EnvMid, this.grid.fluence.tr2EnvMid);

            this.data.vaos.naive.phiEnvLog  = log(abs(this.data.vaos.naive.phiEnvReconNorm));
            this.data.vaos.sp.phiEnvLog  = log(abs(this.data.vaos.sp.phiEnvReconNorm));
%             this.data.vaos.had.phiEnvLog = log(abs(this.data.vaos.had.phiEnvReconNorm));
            
%             this.data.vaos.naive.phiEnvRootLog  = log(abs(this.data.vaos.naiveRoot.phiEnvReconNorm));
%             this.data.vaos.sp.phiEnvRootLog  = log(abs(this.data.vaos.spRoot.phiEnvReconNorm));
%             this.data.vaos.had.phiEnvRootLog = log(abs(this.data.vaos.hadRoot.phiEnvReconNorm));

            if this.oprtn.simSpeckle
                this.data.vaos.speckle.recon.phiSPLog      = log(abs(this.data.vaos.speckle.recon.phiReconSPNorm));
                this.data.vaos.speckle.recon.phiHadLog     = log(abs(this.data.vaos.speckle.recon.phiReconHadNorm));
%                 this.data.vaos.speckleRoot.recon.phiSPLog  = log(abs(this.data.vaos.speckleRoot.recon.phiReconSPNorm));
%                 this.data.vaos.speckleRoot.recon.phiHadLog = log(abs(this.data.vaos.speckleRoot.recon.phiReconHadNorm));
            end
        end

        function generateMathematicalFluence(this)
            this.vaos.createMathematicalFluence();
            %TODO: handle the data of mathematical fluence
        end

        %% Sim functions
        function calcSrcDetFluence(this, i)
            phiSimVars = this.vars.fluence.uVars;
            phiSimVars.mua = this.vars.op.muaVec(i);
            phiSimVars.mus = this.vars.op.musVec(i);
            
            switch this.vars.fluence.LightGeometry
                case 'Point'
                    phiSimVars.mode = 'pointSrc';
                    phiSimVars.srcPos = [0,this.vars.srcParam(1) ,0];
                    
                    this.fs.setVars(phiSimVars)
                    this.fs.config();
                    
                    [this.data.phiSrc, this.data.srcVars] = this.fs.simulate();
                    
                    this.data.phiDet = this.data.phiSrc;
                    this.data.detVars = this.data.srcVars;
                    
                case 'TwoFibers'
                    phiSimVars.mode = 'pointSrc';
                    
                    % Source Fiber:
                    phiSimVars.srcPos = [0,this.vars.fluence.srcParam(1,1) ,0];
                    phiSimVars.srcSize = this.vars.fluence.srcParam(1,2);
                    
                    this.fs.setVars(phiSimVars)
                    this.fs.config();
                    [this.data.phiSrc, this.data.srcVars] = this.fs.simulate();
                    this.savePhi("Src", i);
                    
                    % Detector Fiber:
                    phiSimVars.srcPos = [0,this.vars.fluence.srcParam(2,1) ,0];
                    phiSimVars.srcSize = this.vars.fluence.srcParam(2,2);
                    
                    this.fs.setVars(phiSimVars)
                    this.fs.config();
                    [this.data.phiDet, this.data.detVars] = this.fs.simulate();
                    this.savePhi("Det", i);
                    
                case 'Uniform'
                    phiSimVars.mode = 'Uniform';
                    phiSimVars.srcPos = [0,0,0];
                    
                case 'UniformInf'
                    phiSimVars.mode = 'Uniform';
                    phiSimVars.srcPos = [0,0,0];
                    
                case 'Measured'
                    phiSimVars.mode = 'Measured';
                    phiSimVars.srcPos = [0,0,0];
                    
                    %Src:
                    phiSimVars.pattern = this.data.optics.srcPattern;
                    phiSimVars.patternVec1 = this.data.optics.srcPatternVec1;
                    phiSimVars.patternVec2 = this.data.optics.srcPatternVec2;
                    
                    this.fs.setVars(phiSimVars)
                    this.fs.config();
                    [this.data.fluence.phiSrc{i}, this.vars.fluence.srcVars] = this.fs.simulate();
                    
                    %Det:
                    phiSimVars.pattern = this.data.optics.detPattern;
                    phiSimVars.patternVec1 = this.data.optics.detPatternVec1;
                    phiSimVars.patternVec2 = this.data.optics.detPatternVec2;
                    
                    this.fs.setVars(phiSimVars)
                    this.fs.config();
                    [this.data.fluence.phiDet{i}, this.vars.fluence.detVars] = this.fs.simulate();


                    this.data.fluence.phiSrcNorm{i} =  normMatf(this.data.fluence.phiSrc{i});
                    this.data.fluence.phiDetNorm{i} =  normMatf(this.data.fluence.phiDet{i});
            end

            reset(gpuDevice(1));
        end

        function simFluenceSet(this, uVars)
            this.setSaveVars(uVars)
            this.setFluenceVars(uVars)

            for i=1:this.vars.numOfPhantoms
                this.vars.curSimIdx = i;
               
                fprintf("AOSim: Simulating: mua no. %d \n", i);
                this.calcSrcDetFluence(i);

                this.data.fluence.phiLight{i}     = this.data.fluence.phiSrc{i} .* this.data.fluence.phiDet{i};
                this.data.fluence.phiLightRoot{i} = sqrt(abs(this.data.fluence.phiLight{i}));

                this.data.fluence.phiLightNorm{i} =  this.data.fluence.phiSrcNorm{i} .* this.data.fluence.phiDetNorm{i};
                this.data.fluence.phiLightNormRoot{i} = sqrt(abs(this.data.fluence.phiLightNorm{i}));
            end
            
            if this.vars.save.flag
                this.saveFluence()
            end
            
            this.extractFluenceSim();
            this.extractMuEffFluence();
            this.displayFluenceResults();

            reset(gpuDevice(1));
        end

        function res = simAndAnalyse(this)
%             this.extractFluenceSim();

            if ~this.oprtn.loadVaos
                res.fluence.phi     = this.vaos.interpAlignReplPhi(this.data.fluence.phiEnv, this.grid.fluence.depthVec, true, false, false);
%                 res.fluence.phiRoot = this.vaos.interpAlignReplPhi(this.data.fluence.phiEnvRoot, this.grid.fluence.depthVec, true, true, false);

                fprintf("AOSim: Simulating VAOS on Phi\n");
                res.naive = this.vaos.reconNaive(res.fluence.phi, false);
                res.sp    = this.vaos.reconSP(res.fluence.phi, false);
%                 res.had   = this.vaos.reconHad(res.fluence.phi, false);
                
                if this.oprtn.simSpeckle
                    res.speckle = this.vaos.speckleSim(res.fluence.phi, false);
                end
                
%                 fprintf("AOSim: Simulating VAOS on Rooted Phi\n");
%                 res.naiveRoot = this.vaos.reconNaive(res.fluence.phiRoot, false);
%                 res.spRoot    = this.vaos.reconSP(res.fluence.phiRoot, false);
%                 res.hadRoot   = this.vaos.reconHad(res.fluence.phiRoot, false);

                if this.oprtn.simSpeckle
                    res.speckleRoot = this.vaos.speckleSim(res.fluence.phiRoot, false);
                end
                
                this.data.vaos = res;
                
                this.analyseResults();
                
                if this.vars.save.flag
                    this.vaos.saveVAOSRes(this.vars.save.path)
                end
            end

            this.compareMuEff();
            this.displayResults();

            res.data.meas    = this.data.meas;
            res.data.fluence = this.data.fluence;
            res.data.vaos    = this.data.vaos;

            res.grid.meas    = this.grid.meas;
            res.grid.fluence = this.grid.fluence;
            res.grid.vaos    = this.grid.vaos;
        end
        
        %% Plot and Save
        function displayResults(this)
            % Plot Fluence Sim
            this.displayFluenceResults();

            % Plot Measurement
            isMeas = this.oprtn.loadMeas;
            if isMeas; this.ml.displayResults(); end

            % Plot VAOS Simulations for comparison:
            this.vaos.displayAllRecon(this.data.vaos);

            % Plot All Together (fluence, selectedVAOS, meas):            
            this.displayMuEffComparison();
        end
        
        function displayPulse(this)
            numOfPhantoms = this.vars.numOfPhantoms;
            cols = ceil(sqrt(numOfPhantoms));
            rows = floor(sqrt(numOfPhantoms));
            
            % Plot All Together (fluence, selectedVAOS, meas):            
            phiLog   = this.data.vaos.fluence.phiMidLog;
            naiveLog = this.data.vaos.naive.phiEnvLog;
            convLog  = this.data.vaos.sp.phiEnvLog;

            depthVec        = this.grid.vaos.depthVec; 
            
            legStr = ["Fluence (MCX)", "Naive", "Conv"];

            phSize = this.grid.fluence.depthVec(end);

            xlimLow  = -phSize;
            xlimHigh =  0.1*phSize;
            
            [~, I] = min(abs(depthVec + phSize));

            ylimLow  = min(phiLog(:,I));
            ylimHigh = 0.1;

            isSpec = this.oprtn.simSpeckle;
            if isSpec
                speckleLog = this.data.vaos.speckle.recon.phiSPLog;
                depthVecSpeckle   = this.data.vaos.speckle.vars.xUS;
                legStr = [legStr, "Speckle"];
            end
            
            isMeas = this.oprtn.loadMeas;
            if isMeas
                measLog         = this.data.meas.phiLog;
                depthVecMeas    = this.grid.meas.depthVec;
                legStr = [legStr, "Meas", 'Fit'];
                [~, I] = min(abs(depthVecMeas + phSize));
                ylimLow  = min(measLog(:,I));
            end

            depthVec        = this.grid.vaos.depthVec; 
            
            pulseEnvLog = flip(log(abs(this.data.us.data.pulses.focalPulseEnvNorm)));

            [~,I] = max(pulseEnvLog);
            depthVecPulse = -this.vars.us.grid.pulses.axialVec;
            depthVecPulse = depthVecPulse - depthVecPulse(I);

            figure();
            for i=1:numOfPhantoms
                ax(i) = subplot(rows, cols, i);
                hold(ax(i), 'on')
                plot(ax(i), depthVec, phiLog(i,:))
                plot(ax(i), depthVec, naiveLog(i,:))
                plot(ax(i), depthVec, convLog(i,:))
                if isSpec; plot(ax(i), depthVecSpeckle, speckleLog(i,:)); end
                if isMeas; plot(ax(i), depthVecMeas, measLog(i,:)); end
                plot(ax(i), depthVecPulse, pulseEnvLog)
                xlabel(ax(i), "Depth [mm]");
                ylabel(ax(i), "Fluence [AU - Log]");
                legend(ax(i), legStr, 'Location', 'northwest');
                title(ax(i), sprintf("Phantom-%d", i));
                xlim(ax(i), [xlimLow, xlimHigh])
                ylim(ax(i), [ylimLow ,ylimHigh])
            end
            linkaxes(ax);          
        end
        
        function displayMuEffComparison(this)
            numOfPhantoms = this.vars.numOfPhantoms;
            cols = ceil(sqrt(numOfPhantoms));
            rows = floor(sqrt(numOfPhantoms));
            
            % Plot All Together (fluence, selectedVAOS, meas):            
            phiLog   = this.data.vaos.fluence.phiMidLog;
            naiveLog = this.data.vaos.naive.phiEnvLog;
            convLog  = this.data.vaos.sp.phiEnvLog;

            depthVec        = this.grid.vaos.depthVec; 
            
%             legStr = ["Fluence (MCX)", "Naive", "Conv"];
            legStr = ["Fluence (MCX)", "Conv"];
            phSize = this.grid.fluence.depthVec(end);

            xlimLow  = -phSize;
            xlimHigh =  0.1*phSize;
            
            [~, I] = min(abs(depthVec + phSize));

            ylimLow  = min(phiLog(:,I));
            ylimHigh = 0.1;
            
            muEffData = this.data.muEff.vaos;
            ratio = muEffData.ratio;
            fit = muEffData.fit;
            muEffGT = this.vars.op.muEffGT;

            isSpec = this.oprtn.simSpeckle;
            if isSpec
                speckleLog = this.data.vaos.speckle.recon.phiSPLog;
                depthVecSpeckle   = this.data.vaos.speckle.vars.xUS;
                legStr = [legStr, "Speckle"];
            end
            
            isMeas = this.oprtn.loadMeas;
            if isMeas
                measLog         = this.data.meas.phiLog;
                depthVecMeas    = this.grid.meas.depthVec;
                fitModel        = muEffData.model.meas; 
                for i=1:numOfPhantoms
                    depthVecFit{i,:}     = muEffData.data(i).meas.depthVecCut;
                end
                rsqr            = muEffData.rsqr.meas;
                legStr = [legStr, "Meas", 'Fit'];
                [~, I] = min(abs(depthVecMeas + phSize));
                ylimLow  = min(measLog(:,I));
            end

            figure();
            for i=1:numOfPhantoms
                ax(i) = subplot(rows, cols, i);
                hold(ax(i), 'on')
                
                plot(ax(i), depthVec, phiLog(i,:), 'Color', '#0072BD')
%                 plot(ax(i), depthVec, naiveLog(i,:), 'Color', '#A2142F')
                plot(ax(i), depthVec, convLog(i,:), 'Color', '#D95319')
                if isSpec; plot(ax(i), depthVecSpeckle, speckleLog(i,:), 'Color', '#4DBEEE'); end
                if isMeas
                    plot(ax(i), depthVecMeas, measLog(i,:), 'Color', '#7E2F8E'); 
                    plot(ax(i), depthVecFit{i,:}, fitModel{i,:}, 'Color', '#77AC30')
                end
                xlabel(ax(i), "Depth [mm]");
                ylabel(ax(i), "Fluence [AU - Log]");
                legend(ax(i), legStr, 'Location', 'northwest');
                if isMeas
                    title(ax(i), sprintf("Phantom-%d, R^2=%.3f", i, rsqr(i)));
                else
                    title(ax(i), sprintf("Phantom-%d", i));
                end
                xlim(ax(i), [xlimLow, xlimHigh])
                ylim(ax(i), [ylimLow ,ylimHigh])
            end
            linkaxes(ax);
            
            
            legStr = ["GT", "Fluence", "Conv"];
            if isSpec; legStr = [legStr, "Speckle"]; end
            if isMeas; legStr = [legStr, "Meas"];    end
            figure();
            ax = subplot(1,2,1);
            plot(ax, ratio.gt, ratio.gt, '-o'); hold on;
            plot(ax, ratio.gt, ratio.fluence, '-x'); hold on;
%             plot(ax, ratio.gt, ratio.naive)
            plot(ax, ratio.gt, ratio.conv, '-*' )
            if isSpec;  plot(ax, ratio.gt, ratio.speckle, '-v'); end
            if isMeas;  plot(ax, ratio.gt, ratio.meas, '-^'); end
            xlabel(ax, "GT Ratio");
            ylabel(ax, "Ratio");
            legend(ax, legStr, 'Location', 'northwest');
            title(ax, "In-Set Coeffs. Ratio");

            ax = subplot(1,2,2);
            plot(ax, muEffGT, muEffGT, '-o'); hold on;
            plot(ax, muEffGT, fit.fluence, '-x');
%             plot(ax, muEffGT, fit.naive)
            plot(ax, muEffGT, fit.conv, '-*')
            if isSpec;  plot(ax, muEffGT, fit.speckle, '-v'); end
            if isMeas;  plot(ax, muEffGT, fit.meas, '-^'); end
            xlabel(ax, "GT \mu_{Eff} [mm^{-1}]");
            ylabel(ax, "\mu_{Eff} [mm^{-1}]");
            legend(ax, legStr, 'Location', 'northwest');
            title(ax, "In-Set Coeffs.");
        end
        
        function displayFluenceResults(this)
            numOfPhantoms = this.vars.numOfPhantoms;
            cols = ceil(sqrt(numOfPhantoms));
            rows = floor(sqrt(numOfPhantoms));

            % Plot Fluence Sim
            tr1Mid = this.grid.fluence.tr1Mid;
            tr2Mid = this.grid.fluence.tr2Mid;
            
            depthVec = this.grid.fluence.depthVec;
            
            figure();
            for i=1:numOfPhantoms
                ax = subplot(rows, cols, i);
                hold(ax, 'on')
                plot(ax, depthVec, log(this.data.fluence.phiSrcNorm{i}(:, tr1Mid, tr2Mid)));
                plot(ax, depthVec, log(this.data.fluence.phiDetNorm{i}(:, tr1Mid, tr2Mid)));
                plot(ax, depthVec, log(this.data.fluence.phiLightNormRoot{i}(:, tr1Mid, tr2Mid)));
                xlabel("Depth [mm]");
                ylabel("Fluence [AU - Log]");
                legend("Src", "Det", "Light Root");
                title(sprintf("Phantom-%d", i));
            end
            

%             figure();
%             for i=1:numOfPhantoms
%                 ax = subplot(rows, cols, i);
%                 hold(ax, 'on')
%                 plot(ax, depthVec, this.data.fluence.phiSrc{i}(:, tr1Mid, tr2Mid));
%                 plot(ax, depthVec, this.data.fluence.phiDet{i}(:, tr1Mid, tr2Mid));
%                 plot(ax, depthVec, this.data.fluence.phiLightRoot{i}(:, tr1Mid, tr2Mid));
%                 xlabel("Depth [mm]");
%                 ylabel("Fluence [AU]");
%                 legend("Src", "Det", "Light Root");
%                 title(sprintf("Phantom-%d", i));
%             end
%             
%             figure();
%             for i=1:numOfPhantoms
%                 ax = subplot(rows, cols, i);
%                 hold(ax, 'on')
%                 plot(ax, depthVec, this.data.fluence.phiSrcNorm{i}(:, tr1Mid, tr2Mid));
%                 plot(ax, depthVec, this.data.fluence.phiDetNorm{i}(:, tr1Mid, tr2Mid));
%                 plot(ax, depthVec, this.data.fluence.phiLightNormRoot{i}(:, tr1Mid, tr2Mid));
%                 xlabel("Depth [mm]");
%                 ylabel("Fluence [Normalized AU]");
%                 legend("Src", "Det", "Light Root");
%                 title(sprintf("Phantom-%d", i));
%             end
        end

        function displayResultsForIndexing(this)
            numOfPhantoms = this.vars.numOfPhantoms;
            isMeas = this.oprtn.loadMeas;
            isSpec = this.oprtn.simSpeckle;

            cols = 1 + isMeas + isSpec;
            j = 1;
            figure();
            ax = subplot(1,cols, 1);
            
            % Fluence & Sims:
            phiLog   = this.data.vaos.fluence.phiMidLog;    
            hold(ax, 'on');
            for i=1:numOfPhantoms
                plot(phiLog(i,:));
                xlabel("Depth [Idx]");
                ylabel("Fluence [AU - Log]");
            end
            title("Fluence");

            % Speckle Sim:
            if isSpec
                j = j+1;
                ax = subplot(1, cols, j);
                hold(ax, 'on');
                speckleLog = this.data.vaos.speckle.recon.phiSPLog;
                for i=1:numOfPhantoms
                    plot(speckleLog(i,:));
                    xlabel("Depth [Idx]");
                    ylabel("Fluence [AU - Log]");
                end
                title("Speckle");
            end
                
            % Measurements:
            if isMeas
                j = j+1;
                ax = subplot(1, cols, j);
                hold(ax, 'on');
                measLog = this.data.meas.int.phiLog;
                for i=1:numOfPhantoms
                    plot(measLog(i,:));
                    xlabel("Depth [Idx]");
                    ylabel("Fluence [AU - Log]");
                end
                title("Measurement");
            end
        end

        function saveFluence(this)
            res.vars               = this.vars.fluence;
            res.vars.op            = this.vars.op;
            res.vars.numOfPhantoms = this.vars.numOfPhantoms;
            res.data = this.data.fluence;
            timeStamp = datetime('now','TimeZone','local','Format','d-MMM-y HH-mm-s');
            filename = sprintf("%s/%s-%s.mat", this.vars.save.path, timeStamp, this.vars.save.simName);
            save(filename, '-Struct', 'res', '-v7.3');
        end


    end
end