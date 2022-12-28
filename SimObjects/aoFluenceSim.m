classdef aoFluenceSim <handle
    %AOFLUENCESIM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        usa
        ia
        fs
        cnv
        data
        vars
        paths
        figs
    end
    
    methods (Static)
        function uVars = createUserVars()
            uVars.usPath  = [];
            uVars.srcPath = [];
            uVars.detPath = [];
            
            uVars.phiPath  = [];
            uVars.loadPhi  = false;
            
            uVars.numOfPhantoms = [];
            uVars.geometry      = [];
            uVars.srcParam      = []; %[size, location]
            
            uVars.simulator     = [];

            uVars.muaVec = [];
            uVars.mus    = [];
            uVars.g      = [];
            uVars.ref    = [];
            
            uVars.nphoton = [];
            
            uVars.meshSize = [];
            uVars.meshRes  = [];
            
            uVars.usAx = [];
            uVars.usDir = []; % [1/-1]
            uVars.usFocalDist = []; %distance of the US focus from illumination plane
            uVars.usDataType  = [];
            
            uVars.savePath = '.';
            uVars.saveFlag = false;
            uVars.simName = '';
        end
    end
    
    methods
        function this = aoFluenceSim()
            this.usa = usAnalysis();
            this.fs  = fluenceSim();
            this.ia  = illuminationAnalysis();
            this.cnv = convUsPhi();

            this.paths.curSrcPath = '';
            this.paths.curDetPath = '';
            this.paths.curUSPath = '';
            
            this.vars.table.varNames = ["Sample No. [#]"; "Sim. MuEff [1/mm]"; "Calc. Gradient"; "Calc. MuEff[1/mm]"];
        end
        
        function setVars(this, uVars)
            this.paths.usPath  = uVars.usPath;
            this.paths.src     = uVars.srcPath;
            this.paths.det     = uVars.detPath;
            
            this.paths.phiPath = uVars.phiPath;
            this.vars.loadPhi  = uVars.loadPhi;
            
            % US Analysis:
            this.vars.usDataType = uVars.usDataType;
            
            % Convolution Parameters:
            this.vars.usAx        = uVars.usAx;
            this.vars.usDir       = uVars.usDir;
            this.vars.usFocalDist = uVars.usFocalDist;
            
            % Phi Simulation Parameters:
            this.vars.numOfPhantoms = uVars.numOfPhantom;
            
            this.vars.geometry = uVars.geometry;
            this.vars.srcParam = uVars.srcParam;
            
            this.vars.save.path = uVars.savePath;
            this.vars.save.flag = uVars.saveFlag;
            this.vars.save.simName = uVars.simName;
            
            phiSimVars = fluenceSim.createUserVars();
            
            phiSimVars.simulator = uVars.simulator; % 'mcx', 'toast'
            
            this.vars.muaVec = uVars.muaVec;
            phiSimVars.mua   =  uVars.muaVec(1);
            phiSimVars.mus   = uVars.mus;
            phiSimVars.g     = uVars.g;
            phiSimVars.ref   = uVars.ref;
            
            phiSimVars.nphoton = uVars.nphoton;
            
            phiSimVars.meshSize  = uVars.meshSize; % lenX, lenY, lenZ
            phiSimVars.meshRes   = uVars.meshRes; % [mm]
            
            %Static Configurations:
            phiSimVars.plane       = 'YZ';
            phiSimVars.srcDir      = [1,0,0];
            phiSimVars.loadPattern = false; 
            phiSimVars.intFactor   = 1;
            
            this.vars.phiSimVars = phiSimVars;
            
            % Retreive Simulation Variables:
            fprintf("AOSim: Setting vars to fluence Sim\n");
            this.fs.setVars(phiSimVars)
            this.fs.config();
            this.vars.fluence = this.fs.getVars();

            this.loadDB();
            
            % Set Vars for convolution:
            uVarsCnv = convUsPhi.createUserVars();
            
            uVarsCnv.usAx = uVars.usAx;
            uVarsCnv.usDir = uVars.usDir;
            uVarsCnv.usFocalDist = uVars.usFocalDist;
            uVarsCnv.xVecGrid = this.vars.fluence.xVec;
            uVarsCnv.yVecGrid = this.vars.fluence.yVec;
            uVarsCnv.zVecGrid = this.vars.fluence.zVec;
            
            fprintf("AOSim: Setting vars to Convolution Sim\n");
            this.cnv.setVars(uVarsCnv);
            
            this.vars.cnv = this.cnv.getVars();
            
            this.calcVarsForPP()
        end
        
        function vars = getVars(this)
           vars = this.vars;
        end
        
        function loadDB(this)               
            % Illumination and Detection
            if ~strcmp(this.paths.curSrcPath, this.paths.src)
                fprintf("AOSim: Loading Source Pattern\n");
                [srcPattern, srcVars]=this.ia.analyse(this.paths.src);
                this.data.srcPattern = srcPattern.mask;
                this.data.srcPatternVec1 = srcVars.scanVec-min(srcVars.scanVec);
                this.data.srcPatternVec2 = srcVars.disc1Vec-min(srcVars.disc1Vec);
                this.paths.curSrcPath = this.paths.src;
            end
            
            if ~strcmp(this.paths.curDetPath, this.paths.det)
                fprintf("AOSim: Loading Detector Pattern\n");
                [detPattern, detVars]=this.ia.analyse(this.paths.det);
                this.data.detPattern = detPattern.mask;
                this.data.detPatternVec1 = detVars.scanVec-min(detVars.scanVec);
                this.data.detPatternVec2 = detVars.disc1Vec-min(detVars.disc1Vec);
                this.paths.curDetPath = this.paths.det;
            end
            
            this.ia.flush();
            
            % Ultrasound
            fprintf("AOSim: Loading US Data\n");
            uVarsUS = usAnalysis.createUserVars();
            uVarsUS.usDataPath = this.paths.usPath;
            uVarsUS.usDataType = this.vars.usDataType;
            uVarsUS.intFactor = [1,1,1];
            this.usa.setVars(uVarsUS)
            [this.data.usData,  this.vars.us] = this.usa.analyse();
            this.cnv.setUSData(this.data.usData, this.vars.us);
            this.paths.curUSPath = this.paths.usPath;
        end
        
        function calcVarsForPP(this)
            this.vars.pp.endX = length(this.vars.fluence.xVec);
            this.vars.pp.endY = length(this.vars.fluence.yVec);
            this.vars.pp.endZ = length(this.vars.fluence.zVec);

            this.vars.pp.midX = ceil(this.vars.pp.endX /2);
            this.vars.pp.midY = ceil(this.vars.pp.endX /2);
            this.vars.pp.midZ = ceil(this.vars.pp.endX /2);
            
            this.vars.pp.midYConv = ceil(this.vars.cnv.tr1Size/2);
            this.vars.pp.midZConv = ceil(this.vars.cnv.tr2Size/2);
            
            this.vars.pp.xVecCnv = this.vars.cnv.axVecConv;
            this.vars.pp.dAx = this.vars.us.dAx;
            
            this.vars.pp.startMuEffIdx = find(this.vars.pp.xVecCnv == 0) + 5;
            this.vars.pp.endMuEffIdx   = this.vars.pp.startMuEffIdx + 20;
            
            % Accordin to: "Biomedical Optics Principle and Imaging",
            % L.V.Wang, 2008, pages 98-99 (108-109).
            
            this.vars.pp.simMuEff = sqrt(3*this.vars.muaVec .* (this.vars.muaVec + this.vars.phiSimVars.mus*(1-this.vars.phiSimVars.g) ) );

        end
        
        function postProcessing(this)
            i = this.vars.curSimIdx;
            this.data.midSrc(i,:)   = log(normMatf(this.data.phiSrc{i}(:, this.vars.pp.midY, this.vars.pp.midZ)));
            this.data.midDet(i,:)   = log(normMatf(this.data.phiDet{i}(:, this.vars.pp.midY, this.vars.pp.midZ)));
            this.data.midLight(i,:) = log(normMatf(this.data.phiLight{i}(:, this.vars.pp.midY, this.vars.pp.midZ)));
            
            this.data.midConv(i,:)  = log(normMatf(squeeze(this.data.phiAO{i}(this.vars.pp.midYConv, this.vars.pp.midZConv, :))));
            
            startIdx  = this.vars.pp.startMuEffIdx;
            endIdx = this.vars.pp.endMuEffIdx;
            
%             this.data.finalMuEff(i) = (this.data.midConv(i,startIdx) - this.data.midConv(i,endIdx))/ (this.vars.pp.xVecCnv(startIdx) - this.vars.pp.xVecCnv(endIdx));
            this.data.grad(i) = abs(mean(gradient(this.data.midConv(i,startIdx:endIdx))));
            this.data.calcMuEff(i) = this.data.grad(i) / this.vars.pp.dAx /2;
        end

        function loadPhi(this)
            
        end
        
        function calcFluence(this, i)
            phiSimVars = this.vars.phiSimVars;
            phiSimVars.mua = this.vars.muaVec(i);
            
            switch this.vars.geometry
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
                    phiSimVars.srcPos = [0,this.vars.srcParam(1,1) ,0];
                    phiSimVars.srcSize = this.vars.srcParam(1,2);
                    
                    this.fs.setVars(phiSimVars)
                    this.fs.config();
                    [this.data.phiSrc, this.data.srcVars] = this.fs.simulate();
                    this.savePhi("Src", i);
                    
                    % Detector Fiber:
                    phiSimVars.srcPos = [0,this.vars.srcParam(2,1) ,0];
                    phiSimVars.srcSize = this.vars.srcParam(2,2);
                    
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
                    phiSimVars.pattern = this.data.srcPattern;
                    phiSimVars.patternVec1 = this.data.srcPatternVec1;
                    phiSimVars.patternVec2 = this.data.srcPatternVec2;
                    
                    this.fs.setVars(phiSimVars)
                    this.fs.config();
                    [this.data.phiSrc{i}, this.data.srcVars] = this.fs.simulate();
%                     this.savePhi("Src", i);
                    
                    %Det:
                    phiSimVars.pattern = this.data.detPattern;
                    phiSimVars.patternVec1 = this.data.detPatternVec1;
                    phiSimVars.patternVec2 = this.data.detPatternVec2;
                    
                    this.fs.setVars(phiSimVars)
                    this.fs.config();
                    [this.data.phiDet{i}, this.data.detVars] = this.fs.simulate();
%                     this.savePhi("Det", i);
            end          
        end
        
        function simAndAnalyse(this)
            this.initPlots();
            for i=1:this.vars.numOfPhantoms
                this.vars.curSimIdx = i;
                if this.vars.loadPhi
                    
                    this.loadPhi(i);
                else
                    fprintf("AOSim: Simulating: mua no. %d \n", i);
                    this.calcFluence(i);
                end
                
                this.data.phiLight{i} = this.data.phiSrc{i}.*this.data.phiDet{i};
                fprintf("AOSim: Convolving \n");
                this.data.phiAO{i} = this.cnv.calcConv(this.data.phiLight{i}, this.data.srcVars);
                this.postProcessing()
                this.plotCurrentSim()
            end
           
            this.vars.table.T = table((1:this.vars.numOfPhantoms)', this.vars.pp.simMuEff', this.data.grad', this.data.calcMuEff', 'VariableNames', this.vars.table.varNames);
            display(this.vars.table.T);
            
            if this.vars.save.flag
                this.saveAOSim()
            end
        end
        
        function initPlots(this)
            cols = 3;
            rows = ceil(length(this.vars.muaVec)/cols);
            
            this.figs.hFig = figure();
            
            for i = 1:length(this.vars.muaVec)
                this.figs.ax(i) = subplot(rows, cols, i);
            end
            
            
        end
        
        function plotCurrentSim(this)
            switch this.vars.geometry
                case 'TwoFibers'
                case {'Measured', 'Uniform', 'UniformInf', 'Point'}
                    i = this.vars.curSimIdx;
                    hold(this.figs.ax(i), 'on');
                    this.figs.plt(i,1) = plot(this.figs.ax(i), this.data.srcVars.xVec, this.data.midSrc(i,:));
                    this.figs.plt(i,2) = plot(this.figs.ax(i), this.data.srcVars.xVec, this.data.midDet(i,:));
                    this.figs.plt(i,3) = plot(this.figs.ax(i), this.data.srcVars.xVec, this.data.midLight(i,:));
                    this.figs.plt(i,4) = plot(this.figs.ax(i), this.vars.pp.xVecCnv,   this.data.midConv(i,:));
                    hX = xlabel(this.figs.ax(i), 'X[mm]');
                    hY = ylabel(this.figs.ax(i), 'Phi[mm]');
                    hLeg = legend(this.figs.ax(i), 'Src', 'Det', 'Light', 'Conv'); 
                    set(hLeg, 'Location', 'southwest');
                    hold(this.figs.ax(i), 'off');
            end
        end
        
        function saveAOSim(this)
            simVars = this.getVars();
            simData = this.data;
            simData.usData = simData.usData.pulses;
            timeStamp = strrep(datestr(datetime('now')),':','-');
            filename = sprintf("%s/%s-%s.mat", this.vars.save.path, timeStamp, this.vars.save.simName);
            save(filename, 'simData', 'simVars', '-v7.3');
        end
    end
end

