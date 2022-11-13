classdef toastppSim< handle
    % toastppSim Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %User Parameters
        op
        grid
        meshGrids
        meshGridsHR
        qVec
        
        %ToastObjects
        mesh
        basis
        src
        
        vars
    end
    
    methods (Static)
        function uVars = createUserVars()
            %Mesh & Grid:
            uVars.meshWidth   = []; %mm - half
            uVars.meshDepth   = []; %mm - full
            
            %Optical Properties:
            uVars.ref         = [];
            uVars.musP        = [];
            uVars.muaPivot    = [];
            
            % File System:
            uVars.resDirName  = [];
            uVars.saveFlag = [];
        end
    end
    
    methods
        function this = toastppSim()
            toastThreadCount(20);
            this.vars.curLoadedMesh = '';
        end
        
        function setVars(this, uVars)
            this.op.ref     = uVars.ref;
            this.op.musP    = uVars.musP;
            this.op.mua     = uVars.muaPivot;
            
            this.vars.meshWidth = uVars.meshWidth;
            this.vars.meshDepth = uVars.meshDepth;
        end
        
        function vars = getVars(this)
           vars.grid = this.grid;
           vars.op   = this.op;
           vars.src  = this.src;
           vars.vars = this.vars;
        end
        
        function geometry(this)
            % -------------------------------------------------------------%
            fprintf("Meshing\n");
            this.vars.meshName = sprintf("1LayerSlabMesh-%dx%d-DxW-mm",  this.vars.meshDepth, this.vars.meshWidth);
            meshPath = "./Meshes";
            meshFilename = sprintf("%s.msh", this.vars.meshName);
            meshFullName = sprintf("%s/%s", meshPath, meshFilename);

            if strcmp(this.vars.meshName, this.vars.curLoadedMesh)
                fprintf("Mesh Already Loaded\n");
                return
            end
            
            %Create the mesh from the geometrical points defined in the .geo file
            fprintf("Constructing Mesh\n"); 
            if ~exist(meshFullName, 'file')
                fprintf("Can't locate the mesh. Generating...\n"); 
                str(1) = sprintf("depth = %d;", meshDepth);
                str(2) = sprintf("width = %d;", meshWidth);
                replaceLineInFile("./1LayerBox.geo", [1,2],  str);
                system (sprintf("gmsh -3 ./1LayerBox.geo -o %s", meshFullName));
            else
                fprintf("Located the mesh. Loading without generating.\n");
            end

            % Load the mesh into toastMesh object.
            % Notice: If the .msh file was already constructed the above code is
            % redundant and the mesh should only be loaded with the code below.
            fprintf("Loading Mesh\n");
            this.mesh = toastMesh(meshFullName,'gmsh');
            this.vars.curLoadedMesh = this.vars.meshName;
            fprintf("Done Loading Mesh\n");

            % -------------------------------------------------------------%
            fprintf("Gridding\n");
            
            dl  = 0.5; % distance between 2 elements in the mesh in each one of the axes.

            bbMesh = this.mesh.BoundingBox;

            %Define axes lengths:
            xLen    = bbMesh(2,1) - bbMesh(1,1);  % length of x axis of the mesh -[mm]
            yLen    = bbMesh(2,2) - bbMesh(1,2);  % length of y axis of the mesh -[mm] 
            zLen    = bbMesh(2,3) - bbMesh(1,3);  % length of z axis of the mesh -[mm]

            % Direction with relation to the illumination:
            % Lateral axis      - X [0, len]
            % Transversal plane - YZ [-len/2, len/2]

            % Create axes vectors:
            xVec = bbMesh(1,1):dl:bbMesh(2,1);
            yVec = bbMesh(1,2):dl:bbMesh(2,2);
            zVec = bbMesh(1,3):dl:bbMesh(2,3);

            % Calculate Indexes:
            xSize = length(xVec);
            ySize = length(yVec);
            zSize = length(zVec);
            
            grd = [xSize, ySize, zSize];
            
            % Index of the axis center:
            xMidIdx = floor(xSize/2); 
            yMidIdx = floor(ySize/2);
            zMidIdx = floor(zSize/2);

            % Create mesh grid that will help cast the toastMesh onto:
            [this.meshGrids{1},this.meshGrids{2},this.meshGrids{3}]  = meshgrid(yVec, xVec, zVec);

            % Create high resolution grid high resolution simulation:
            dlHR = dl/10;
            xVecHR = xVec;
            yVecHR = bbMesh(1,2):dlHR:bbMesh(2,2);
            zVecHR = bbMesh(1,3):dlHR:bbMesh(2,3);

            xSizeHR = length(xVecHR);
            ySizeHR = length(yVecHR);
            zSizeHR = length(zVecHR); 
            
            grdHR = [xSizeHR, ySizeHR, zSizeHR];
            
            xMidHRIdx = floor(xSizeHR/2)+1;
            yMidHRIdx = floor(ySizeHR/2)+1;
            zMidHRIdx = floor(zSizeHR/2)+1;
            
            [this.meshGridsHR{1},this.meshGridsHR{2},this.meshGridsHR{3}] = meshgrid(yVecHR, xVecHR, zVecHR);
            
%             curVars = who();
%             for i = 1:length(curVars)
%                 if strcmp(curVars{i}, 'this'); continue; end
%                 this.grid.(curVars{i}) = eval(curVars{i});
%             end

            % -------------------------------------------------------------%
            fprintf("Basising\n");
            tic
            this.basis = toastBasis(this.mesh, this.grid.grd);
            toc

            this.vars.ne             = this.mesh.ElementCount;
            this.vars.nv             = this.mesh.NodeCount;
            this.vars.regidx         = this.mesh.Region;
            this.vars.regno          = unique(this.vars.regidx);%sorted from small to largest.

            this.vars.firstLayerIdx  = find(this.vars.regidx == this.vars.regno(1));
        end

        function setOpticalProperties(this, curRef, curMusP, curMua)
            this.curOp.ref = curRef;
            this.curOp.mua = curMusP;
            this.curOp.mus = curMua;
        end
        
        function res = calcPhi(this)
            % Set the current source:
            this.mesh.SetQM(this.curSrc.pos, []);
            qVec = this.mesh.Qvec('Neumann',this.src.srcType, this.curSrc.size);
            
            % Set optical properties to mesh:
            ref = ones(this.vars.ne,1)*this.op.ref;
            mus = this.curOp.mus*ones(this.vars.ne,1);
            mua = this.curOp.mua*ones(this.vars.ne,1);
  
            % Creating The System Matrix
            fprintf("Creating SysMat\n");
            tic
            K = dotSysmat(this.mesh, mua ,mus, ref, 'EL');
            toc

            % Solve
            fprintf("Solving\n");
            tic
            phi = K\qVec;
            phiMapped  = this.basis.Map('M->B',phi);
            toc
            
            fprintf("Analyzing\n");
            % Extract Phi: Linear and Log
            phiGrd = reshape(phiMapped, this.grid.xSize, this.grid.ySize, this.grid.zSize);
            
            fprintf("Interpolating\n");
            tic
%             phiHR = interp3(X, Y, Z, phiGrd, XHR, YHR, ZHR);
            phiHR = interp3(this.grid.yVec, this.grid.xVec', this.grid.zVec,...
                            phiGrd,...
                            this.grid.yVecHR, this.grid.xVecHR', this.grid.zVecHR);
            toc
            
            fprintf("Extracting\n");

            res.phi      = phiGrd;         %3D
            res.phiHR    = phiHR;
        end
        
    end
end



