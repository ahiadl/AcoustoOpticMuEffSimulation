classdef fluenceSim < handle
    % Simulates fluence based on MCX/Toast++ tools.
    % Illumination is always from YZ plane and towards X direction
    
    properties
        ia
        sim
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

            uVars.objParams = [];
        end
    end
    
    methods
        function this = fluenceSim(type)
            this.ia  = illuminationAnalysis();
            this.vars.simulator = type;
            switch type
                case 'MCX'
                    this.sim = mcxSim();
                    
                case 'Taost++'
                    this.sim = toastppSim();
            end
        end
        
        function setVars(this, uVars)
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

            this.vars.objParams = uVars.objParams;
            
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
           uVarsSim = this.sim.createUserVars();

           switch this.vars.simulator
                case 'MCX'
                    this.vars.xLims =  [0,this.vars.meshSize(2)];
                    this.vars.yLims = [-this.vars.meshSize(2)/2,  this.vars.meshSize(2)/2];
                    this.vars.zLims = [-this.vars.meshSize(3)/2,  this.vars.meshSize(3)/2];
                    
                    uVarsSim.xLims = this.vars.xLims;
                    uVarsSim.yLims = this.vars.yLims;
                    uVarsSim.zLims = this.vars.zLims;

                    uVarsSim.res       = this.vars.meshRes;
                    uVarsSim.intFactor = this.vars.intFactor;
                    
                    switch this.vars.mode
                        case 'pointSrc'
                            uVarsSim.srcType = 'pencil';                           
                        case 'Uniform'
                            uVarsSim.srcType = 'Planar';
                            uVarsSim.planePoint1   = [61, 0, 0];
                            uVarsSim.planePoint2   = [0, 61, 0];
                            uVarsSim.plane   = this.vars.plane;
                        case 'Measured'
                            uVarsSim.srcType = 'pattern';
                            uVarsSim.plane   = this.vars.plane;
                            if this.vars.loadPattern
                               [ptrn, ptrnVars] = this.ia.analyse(this.vars.patternPath);
                               this.vars.pattern  = ptrn.mask;
                               this.vars.patternVec1 = ptrnVars.scanVec - min(ptrnVars.scanVec);
                               this.vars.patternVec2 = ptrnVars.disc1Vec - min(ptrnVars.disc1Vec);
                            end
                            
                            uVarsSim.pattern =  this.vars.pattern;
                            uVarsSim.patternVec1 = this.vars.patternVec1;
                            uVarsSim.patternVec2 = this.vars.patternVec2;
                            uVarsSim.plane   = this.vars.plane;
                        case 'other'
                        otherwise
                            error("Invalid source type.\n");
                    end
                    
                    %Static src configurations:
                    uVarsSim.srcPos  = [0,0,0];
                    uVarsSim.srcDir  = [1 0 0]; 

                    % Optical Properties:
                    uVarsSim.mua = this.vars.mua;
                    uVarsSim.mus = this.vars.mus;
                    uVarsSim.g   = this.vars.g;
                    uVarsSim.ref = this.vars.ref;

                    % Static configurations:
                    uVarsSim.numOfObjects = 1;
                    uVarsSim.cw = true;
                    uVarsSim.tStart = 0;
                    uVarsSim.tstep  = 1e-10;
                    uVarsSim.tend   = 5e-9;
            
                    uVarsSim.nphoton = this.vars.nphoton; 
                    
                    uVarsSim.objParams = this.vars.objParams;

                case 'toast'
                    % TODO: Implement Support.
           end
            
            this.sim.setVars(uVarsSim);
                    
            this.vars.varsSim = this.sim.getVars();
            
            this.vars.xVec = this.vars.varsSim.grid.xVec;
            this.vars.yVec = this.vars.varsSim.grid.yVec;
            this.vars.zVec = this.vars.varsSim.grid.zVec;
        end
        
        function vars = getVars(this)
            vars = this.vars;
        end
        
        function [res,vars] = simulate(this)
            vars = this.vars;

            res = this.sim.calcPhi();

%             switch this.vars.simulator
%                 case 'mcx'
%                     
%                 case 'toast'
% 
%             end
        end
    end
end

