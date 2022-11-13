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
            this.vars.fs = this.fs.getVars();

            this.loadDB();
            
            % Set Vars for convolution:
            uVarsCnv = convUsPhi.createUserVars();
            
            uVarsCnv.usAx = uVars.usAx;
            uVarsCnv.usDir = uVars.usDir;
            uVarsCnv.usFocalDist = uVars.usFocalDist;
            uVarsCnv.xVecGrid = this.vars.fs.xVec;
            uVarsCnv.yVecGrid = this.vars.fs.yVec;
            uVarsCnv.zVecGrid = this.vars.fs.zVec;
            
            
            fprintf("AOSim: Setting vars to Convolution Sim\n");
            this.cnv.setVars(uVarsCnv);
            
            this.vars.cnv = this.cnv.getVars();
            
            this.calcVarsForPP()
        end
        
        function calcVarsForPP(this)
            this.vars.pp.endX = length(this.vars.fs.xVec);
            this.vars.pp.endY = length(this.vars.fs.yVec);
            this.vars.pp.endZ = length(this.vars.fs.zVec);

            this.vars.pp.midX = ceil(this.vars.pp.endX /2);
            this.vars.pp.midY = ceil(this.vars.pp.endX /2);
            this.vars.pp.midZ = ceil(this.vars.pp.endX /2);
            
            this.vars.pp.midYConv = ceil(this.vars.cnv.tr1Size/2);
            this.vars.pp.midZConv = ceil(this.vars.cnv.tr2Size/2);
            
            this.vars.pp.xVecCnv = this.vars.cnv.axVecConv;
        end
        
        function postProcessing(this)
            this.data.midSrc = log(normMatf(this.data.phiSrc(:, this.vars.pp.midY, this.vars.pp.midZ)));
            this.data.midDet = log(normMatf(this.data.phiDet(:, this.vars.pp.midY, this.vars.pp.midZ)));
            this.data.midLight = log(normMatf(this.data.phiLight(:, this.vars.pp.midY, this.vars.pp.midZ)));
            
            this.data.midConv = log(normMatf(squeeze(this.data.phiAO(this.vars.pp.midYConv, this.vars.pp.midZConv, :))));
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
%             if ~strcmp(this.paths.curUSPath, this.paths.usPath)
                fprintf("AOSim: Loading US Data\n");
                uVarsUS = usAnalysis.createUserVars();
                uVarsUS.usDataPath = this.paths.usPath;
                uVarsUS.usDataType = this.vars.usDataType;
                uVarsUS.intFactor = [1,1,1];
                this.usa.setVars(uVarsUS)
                [this.data.usData,  this.vars.us] = this.usa.analyse();
                this.cnv.setUSData(this.data.usData, this.vars.us);
                this.paths.curUSPath = this.paths.usPath;
%             end
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
                    [this.data.phiSrc, this.data.srcVars] = this.fs.simulate();
                    this.savePhi("Src", i);
                    
                    %Det:
                    phiSimVars.pattern = this.data.detPattern;
                    phiSimVars.patternVec1 = this.data.detPatternVec1;
                    phiSimVars.patternVec2 = this.data.detPatternVec2;
                    
                    this.fs.setVars(phiSimVars)
                    this.fs.config();
                    [this.data.phiDet, this.data.detVars] = this.fs.simulate();
                    this.savePhi("Det", i);
            end          
        end
        
        function simAndAnalyse(this)
            for i=1:this.vars.numOfPhantoms
                if this.vars.loadPhi
                    
                    this.loadPhi(i);
                else
                    fprintf("AOSim: Simulating: mua no. %d \n", i);
                    this.calcFluence(i);
                end
                
                this.data.phiLight = this.data.phiSrc.*this.data.phiDet;
                fprintf("AOSim: Convolving \n");
                this.data.phiAO = this.cnv.calcConv(this.data.phiLight, this.data.srcVars);
                this.postProcessing()
                this.saveAOSim(i)
                this.plotCurrentSim()
            end
        end
        
        function plotCurrentSim(this)
            switch this.vars.geometry
                case 'TwoFibers'
                case {'Measured', 'Uniform', 'UniformInf', 'Point'}
                    figure()
                    subplot(3,3,1)
                    imagesc(this.data.srcVars.yVec, this.data.srcVars.xVec, log(squeeze(this.data.phiSrc(:,:,this.vars.pp.midZ))))
                    xlabel( 'Y [mm]' )
                    ylabel( 'X [mm]' )
                    axis tight equal
      
                    subplot(3,3,2)
                    imagesc(this.data.detVars.yVec, this.data.detVars.xVec, log(squeeze(this.data.phiDet(:,:,this.vars.pp.midZ))))
                    xlabel( 'Y [mm]' )
                    ylabel( 'X [mm]' )
                    axis tight equal
     
                    subplot(3,3,3)
                    imagesc(this.data.srcVars.yVec, this.data.srcVars.xVec, log(squeeze(this.data.phiLight(:,:,this.vars.pp.midZ))))
                    xlabel( 'Y [mm]' )
                    ylabel( 'X [mm]' )
                    axis tight equal

                    subplot(3,3,7)
                    plot(this.data.srcVars.xVec, this.data.midSrc); hold on
                    plot(this.data.srcVars.xVec, this.data.midDet);
                    plot(this.data.srcVars.xVec, this.data.midLight);
                    plot(this.vars.pp.xVecCnv, this.data.midConv);
                    xlabel('X[mm]')
                    ylabel('Phi[mm]')
                    legend('Src', 'Det', 'Light');
                    
%                     subplot(3,3,4)
%                     hs=slice(log(squeeze(this.data.phiSrc)), [1,midX,endX],[1, midY, endY],[1,midZ endZ]);
%                     set(hs,'linestyle','none');
%                     axis tight equal
%                     subplot(3,3,5)
%                     hs=slice(log(squeeze(this.data.phiDet)), [1,midX,endX],[1, midY, endY],[1,midZ endZ]);
%                     set(hs,'linestyle','none');
%                     axis tight equal
%                     subplot(3,3,6)
%                     hs=slice(log(squeeze(this.data.phiSrc)), [1,midX,endX],[1, midY, endY],[1,midZ endZ]);
%                     set(hs,'linestyle','none');
%                     axis tight equal
                
            end
        end
        
        function savePhi(this, name, idx)
            
        end
        
        function saveAOSim(this, idx)
            
        end
    end
end

