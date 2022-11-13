classdef fluenceSim < handle
    %FLUENCESIM Summary of this class goes here
    %   Detailed explanation goes here
    % illumination is always from YZ plane and towards X direction
    
    properties
        ia
        simMCX
        simToast
        vars
        data
    end
    
    methods (Static)
        function uVars = createUserVars()
            uVars.simulator = 'mcx'; % 'mcx', 'toast'
            uVars.mode      = 'pointSrc'; % twoFibers; Uniform; Measured
            uVars.meshMode  = []; % mcxOnly: 'slab', 'semiInf'
            uVars.meshSize  = []; % lenX, lenY, lenZ
            uVars.srcPos    = [0,0,0]; % [X,Y,Z] in mm; In case of two fibers this should be vector that specifices distance and angle
            uVars.srcSize   = []; % [mm]
            uVars.meshRes   = []; % [mm]
            uVars.intFactor = []; %
            uVars.srcDir    = []; % direction, as a normalized vector
            
            uVars.mua = []; %
            uVars.mus = []; %
            uVars.g   = []; %
            uVars.ref = []; %
            
            uVars.loadPattern = false;
            uVars.patternPath = [];
            uVars.pattern = [];
            uVars.patternVec1 = [];
            uVars.patternVec2 = [];
            
            uVars.nphoton = [];
            
            uVars.plane = [];
        end
    end
    
    methods
        function this = fluenceSim()
            this.ia  = illuminationAnalysis();
            this.simMCX = mcxSim();
            this.simToast = toastppSim();
        end
        
        function setVars(this, uVars)
            this.vars.simulator = uVars.simulator;
            this.vars.mode      = uVars.mode;     % twoFibers; Uniform; Measured
            this.vars.meshMode  = uVars.meshMode; % mcxOnly: 'slab', 'semiInf'
            this.vars.meshSize  = uVars.meshSize; % [lenX, lenY, lenZ]
            this.vars.srcPos    = uVars.srcPos;   % [X,Y,Z] - in case of two fibers this should be vector that specifices distance and angle
            this.vars.srcSize   = uVars.srcSize;
            this.vars.meshRes   = uVars.meshRes;
            this.vars.intFactor = uVars.intFactor;
            
            this.vars.srcDir = uVars.srcDir;   % mcxOnly
            
            this.vars.mua  = uVars.mua;
            this.vars.mus  = uVars.mus;
            this.vars.g    = uVars.g;
            this.vars.musP = this.vars.mus * (1-this.vars.g);
            this.vars.ref  = uVars.ref;
            
            this.vars.loadPattern = uVars.loadPattern;
            
            if strcmp(this.vars.mode, 'Measured')
                if this.vars.loadPattern
                    this.vars.patternPaths = uVars.patternPath;
                else
                    this.vars.pattern = uVars.pattern;
                    this.vars.patternVec1 = uVars.patternVec1;
                    this.vars.patternVec2 = uVars.patternVec2;
                end
            end
            
            if strcmp(this.vars.mode, 'Uniform') || strcmp(this.vars.mode, 'Measured')
               this.vars.plane = uVars.plane; 
            end
            
            
            %MCX Only:
            this.vars.nphoton = uVars.nphoton;
            
        end
        
        function config(this)
           switch this.vars.simulator
                case 'mcx'
                    uVarsMcx = this.simMCX.createUserVars();
                    
                    this.vars.xLims =  [0,this.vars.meshSize(2)];
                    this.vars.yLims = [-this.vars.meshSize(2)/2,  this.vars.meshSize(2)/2];
                    this.vars.zLims = [-this.vars.meshSize(3)/2,  this.vars.meshSize(3)/2];
                    
                    uVarsMcx.xLims = this.vars.xLims;
                    uVarsMcx.yLims = this.vars.yLims;
                    uVarsMcx.zLims = this.vars.zLims;

                    uVarsMcx.res       = this.vars.meshRes;
                    uVarsMcx.intFactor = this.vars.intFactor;
                    
                    switch this.vars.mode
                        case 'pointSrc'
                            uVarsMcx.srcType = 'pencil';                           
                        case 'Uniform'
                            uVarsMcx.srcType = 'Planar';
                            uVarsMcx.planePoint1   = [61, 0, 0];
                            uVarsMcx.planePoint2   = [0, 61, 0];
                            uVarsMcx.plane   = this.vars.plane;
                        case 'Measured'
                            uVarsMcx.srcType = 'pattern';
                            uVarsMcx.plane   = this.vars.plane;
                            if this.vars.loadPattern
                               [ptrn, ptrnVars] = this.ia.analyse(this.vars.patternPath);
                               this.vars.pattern  = ptrn.mask;
                               this.vars.patternVec1 = ptrnVars.scanVec - min(ptrnVars.scanVec);
                               this.vars.patternVec2 = ptrnVars.disc1Vec - min(ptrnVars.disc1Vec);
                            end
                            
                            uVarsMcx.pattern =  this.vars.pattern;
                            uVarsMcx.patternVec1 = this.vars.patternVec1;
                            uVarsMcx.patternVec2 = this.vars.patternVec2;
                            uVarsMcx.plane   = this.vars.plane;
                        case 'other'
                        otherwise
                            error("Invalid source type.\n");
                    end
                    
                    %Static src configurations:
                    uVarsMcx.srcPos  = [0,0,0];
                    uVarsMcx.srcDir  = [1 0 0]; 

                    % Optical Properties:
                    uVarsMcx.mua = this.vars.mua;
                    uVarsMcx.mus = this.vars.mus;
                    uVarsMcx.g   = this.vars.g;
                    uVarsMcx.ref = this.vars.ref;

                    % Static configurations:
                    uVarsMcx.numOfObjects = 1;
                    uVarsMcx.cw = true;
                    uVarsMcx.tStart = 0;
                    uVarsMcx.tstep  = 1e-10;
                    uVarsMcx.tend   = 5e-9;
            
                    uVarsMcx.nphoton = this.vars.nphoton; 
                    
                    this.simMCX.setVars(uVarsMcx);
                    
                    this.vars.varsSim = this.simMCX.getVars();
                    
                    this.vars.xVec = this.vars.varsSim.grid.xVec;
                    this.vars.yVec = this.vars.varsSim.grid.yVec;
                    this.vars.zVec = this.vars.varsSim.grid.zVec;
                case 'toast'
                    
                    
           end
        end
        
        function vars = getVars(this)
            vars = this.vars;
        end
        
        function [res,vars] = simulate(this)
            vars = this.vars;
            switch this.vars.simulator
                case 'mcx'
                    res = this.simMCX.calcPhi();
                case 'toast'

            end
        end
    end
end

