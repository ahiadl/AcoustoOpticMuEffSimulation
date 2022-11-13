classdef mcxSim <handle
    %MCXSIM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cfg
        src
        det
        grid
        mesh
        op
        time
        mcx
    end
    
    methods (Static)
        function uVars = createUserVars()
            uVars.xLims = [];
            uVars.yLims = [];
            uVars.zLims = [];
            
            uVars.res = [];
            
            uVars.srcType = 'pencil';
            uVars.srcPos  = [];
            uVars.srcDir  = [];
            uVars.pattern = [];
            
            uVars.cyclicBoundary = false;
            
            uVars.patternWidth  = [];  %%
            uVars.patternHeight = [];  %%
            uVars.waistSize     = [];  %%
            uVars.halfAngle     = [];  %%
            uVars.lineEndPoint  = [];  %%
            uVars.slitEndPoint  = [];  %%
            
            
            uVars.patternSizeFactor = 1;
            uVars.patternVec1 = [];
            uVars.patternVec2 = [];
            
            uVars.plane = '';
            uVars.planeSize = [];
            
            uVars.tstart = 0;
            uVars.tstep  = 1e-10;
            uVars.tend   = 5e-9;
            
            uVars.cw = true;
            
            uVars.numOfObjects = []; %TODO: maybe should be in geometry?
            
            uVars.mua = [];
            uVars.mus = [];
            uVars.g   = [];
            uVars.ref = [];
            
            uVars.nphoton = 1e9;
        end
    end
    
    methods
        function this = mcxSim()
            this.mcx.autopilot  = 1; %1-automatically set threads and blocks, [0]-use nthread/nblocksize
            this.mcx.nblocksize = 64; %how many CUDA thread blocks to be used [64]
            this.mcx.nthread    = 4096; %the total CUDA thread number [2048]
        end
        
        function setVars(this, uVars)
            % ------ Geometry (Grid & Mesh) ---------
            this.grid.xLims = uVars.xLims;
            this.grid.yLims = uVars.yLims;
            this.grid.zLims = uVars.zLims;
            
            this.grid.res   = uVars.res;
            
            this.grid.xVec =  this.grid.xLims(1):this.grid.res:this.grid.xLims(2);
            this.grid.yVec =  this.grid.yLims(1):this.grid.res:this.grid.yLims(2);
            this.grid.zVec =  this.grid.zLims(1):this.grid.res:this.grid.zLims(2);
            
            this.grid.size = [length(this.grid.xVec), length(this.grid.yVec), length(this.grid.zVec)];
            
            this.grid.cyclicBoundary = uVars.cyclicBoundary;
            
            this.geometry();
            
            % ------ Source ---------
            this.src.type = uVars.srcType;
            
            this.src.userPos = uVars.srcPos;
            this.src.pos  = uVars.srcPos - [this.grid.xVec(1), this.grid.yVec(1), this.grid.zVec(1)]/this.grid.res+1;
            
            this.src.dir  = uVars.srcDir;
            this.src.pattern = uVars.pattern;
            this.src.numPattern = size(this.src.pattern,3);

            this.src.patternVec1 = uVars.patternVec1;
            this.src.patternVec2 = uVars.patternVec2;
            
            this.src.plane = uVars.plane;
            this.src.patternSizeFactor =uVars.patternSizeFactor;
            
            this.src.planeSize     = uVars.planeSize;
            this.src.patternWidth  = uVars.patternWidth;
            this.src.patternHeight = uVars.patternHeight;
            this.src.waistSize     = uVars.waistSize;
            this.src.halfAngle     = uVars.halfAngle;
            this.src.lineEndPoint  = uVars.lineEndPoint;
            this.src.slitEndPoint  = uVars.slitEndPoint;
            
            this.src.nphoton = uVars.nphoton;

            this.setSource();
            
            % ------ Detector ---------

            this.setDetector();
            
            % ------ Optical Properties ---------
            this.op.numOfObjects = uVars.numOfObjects; %TODO: maybe should be in geometry?
            
            % These parameters should be vectors of length uVars.numOfObjects
            this.op.mua = uVars.mua;
            this.op.mus = uVars.mus;
            this.op.g   = uVars.g;
            this.op.ref = uVars.ref;
            
            this.setOpticalProperties();
            
            % ------ Time ---------
            this.time.tstart = uVars.tstart;
            this.time.tstep  = uVars.tstep;
            this.time.tend   = uVars.tend;
            
            this.time.cw     = uVars.cw;
            
            this.setTime();

            % ------ Simulation ---------
            
            this.setSimulation();
        end
        
        function vars = getVars(this)
            vars.grid = this.grid;
            vars.op = this.op;
            vars.src = this.src;
            vars.time = this.time;
        end
        
        function geometry(this)
            %% Parameters:
            % cfg.vol: a 3D array specifying the media index in the domain.
            %           can be uint8, uint16, uint32, single or double arrays.
            %           2D simulations are supported if cfg.vol has a singleton
            %           dimension (in x or y); srcpos/srcdir must belong to
            %           the 2D plane in such case.
            %           MCXLAB also accepts 4D arrays to define continuously varying media. 
            %           The following formats are accepted
            %            1 x Nx x Ny x Nz float32 array: mua values for each voxel (must use permute to make 1st dimension singleton)
            %            2 x Nx x Ny x Nz float32 array: mua/mus values for each voxel (g/n use prop(2,:))
            %            4 x Nx x Ny x Nz uint8   array: mua/mus/g/n gray-scale (0-255) interpolating between prop(2,:) and prop(3,:)
            %            2 x Nx x Ny x Nz uint16  array: mua/mus gray-scale (0-65535) interpolating between prop(2,:) and prop(3,:)
            % cfg.unitinmm: defines the length unit for a grid edge length [1.0]            
            % cfg.bc: per-face boundary condition (BC), a strig of 6 letters (case insensitive) for
            %          bounding box faces at -x,-y,-z,+x,+y,+z axes;
            %          overwrite cfg.isreflect if given.
            %          each letter can be one of the following:
            %          '_': undefined, fallback to cfg.isreflect
            %          'r': like cfg.isreflect=1, Fresnel reflection BC
            %          'a': like cfg.isreflect=0, total absorption BC
            %          'm': mirror or total reflection BC
            %          'c': cyclic BC, enter from opposite face
            %          in addition, cfg.bc can contain up to 12 characters,
            %          with the 7-12 characters indicating bounding box
            %          facets -x,-y,-z,+x,+y,+z are used as a detector. The 
            %          acceptable characters for digits 7-12 include
            %          '0': this face is not used to detector photons
            %          '1': this face is used to capture photons (if output detphoton)
            this.cfg.vol = ones(this.grid.size, 'uint8');
            this.cfg.unitinmm = this.grid.res;
            if this.grid.cyclicBoundary
                 this.cfg.bc  = 'ccrccr'; %TODO: add logic to claculate which boundaries should be chosen
            else
                 this.cfg.bc  = '______';
            end
        end
        
        function setSource(this)
            %% Parameters:
            % cfg.srcpos:     a 1 by 3 vector, the position of the source in grid unit
            % cfg.srcdir:     a 1 by 3 vector, specifying the incident vector; if srcdir
            %                 contains a 4th element, it specifies the focal length of
            %                 the source (only valid for focuable src, such as planar, disk,
            %                 fourier, gaussian, pattern, slit, etc); if the focal length
            %                 is nan, all photons will be launched isotropically regardless
            %                 of the srcdir direction.
            % cfg.issrcfrom0: 1-first voxel is [0 0 0], [0]- first voxel is [1 1 1]

            this.cfg.nphoton = this.src.nphoton;
            
            this.cfg.srcpos = this.src.pos;
            this.cfg.srcdir = this.src.dir;
            this.cfg.issrcfrom0=0;
            
            switch this.src.type
                case 'pencil'
                    % 'pencil' - default, pencil beam, no param needed                    
                    this.cfg.srctype = 'pencil';
                case 'isotropic'
                    % 'isotropic' - isotropic source, no param needed
                    this.cfg.srctype = 'isotropic';
                case 'cone'
                    % 'cone' - uniform cone beam, srcparam1(1) is the 
                    % half-angle in radian
                    this.cfg.srctype = 'cone';
                    this.cfg.srcparam1 = halfAngle;
                case 'gaussian'
                    % 'gaussian' [*] - a collimated gaussian beam, srcparam1(1) 
                    % specifies the waist radius (in voxels)
                    this.cfg.srctype = 'gaussian';
                    this.cfg.srcparam1 = waistSize;
                case 'planar'
                    % 'planar' [*] - a 3D quadrilateral uniform planar 
                    % source, with three corners specified by srcpos, 
                    % srcpos+srcparam1(1:3) and srcpos+srcparam2(1:3)
                    this.cfg.srctype = 'planar';
                    planeSizeIdx = this.src.planeSize/this.grid.res;
%                     this.cfg.srcparam1 = this.src.planePoint1;
%                     this.cfg.srcparam2 = this.src.planePoint2;
                    switch this.src.plane
                        case 'XY'
                            this.cfg.srcpos    =  [this.src.pos(1:2) -  planeSizeIdx/2, 0];  
                            this.cfg.srcparam1 = [planeSizeIdx(1), 0, 0,  size(this.src.pattern,1)];
                            this.cfg.srcparam2 = [0, planeSizeIdx(2), 0,  size(this.src.pattern,2)];
                        case 'YZ'
                            this.cfg.srcpos    = [0, this.src.pos(2:3) - planeSizeIdx/2];  
                            this.cfg.srcparam1 = [0, planeSizeIdx(1), 0,  size(this.src.pattern,1)];
                            this.cfg.srcparam2 = [0, 0, planeSizeIdx(2),  size(this.src.pattern,2)];
                        case 'XZ'
                            this.cfg.srcpos    = [this.src.pos(1)- planeSizeIdx(1)/2, 0, this.src.pos(3) - planeSizeIdx(2)/2];  
                            this.cfg.srcparam1 = [planeSizeIdx(1), 0, 0,  size(this.src.pattern,1)];
                            this.cfg.srcparam2 = [0, 0, planeSizeIdx(2),  size(this.src.pattern,2)];
                        case 'custom'                            
                            this.cfg.srcparam1=[this.src.planePoint1, size(this.src.pattern,1)]; %srcparam1(1:3)=plane corner, srcparam1(4) pattern Dimension X
                            this.cfg.srcparam2=[this.src.planePoint2, size(this.src.pattern,2)]; %srcparam2(1:3)=plane corner, srcparam2(4) pattern Dimension Y
                    end
  
                case 'pattern'
                    % 'pattern' [*] - a 3D quadrilateral pattern illumination, same as above, except
                    % srcparam1(4) and srcparam2(4) specify the pattern array x/y dimensions,
                    % and srcpattern is a floating-point pattern array, with values between [0-1]. 
                    % if cfg.srcnum>1, srcpattern must be a floating-point array with 
                    % a dimension of [srcnum srcparam1(4) srcparam2(4)]
                    % Example: <demo_photon_sharing.m>
                    % cfg.srcpattern: see cfg.srctype for details
                    % cfg.srcnum: the number of source patterns that are
                    %   simultaneously simulated; only works for 'pattern'
                    %   source, see cfg.srctype='pattern' for details
                    %   Example <demo_photon_sharing.m>

                    this.cfg.srctype    = 'pattern';
                    this.cfg.srcpattern = this.src.pattern; %2D pattern
                    this.cfg.srcnum     = 1;
                    patternLen          =  [this.src.patternVec1(end), this.src.patternVec2(end)]...
                                           * this.src.patternSizeFactor;
                    patternGridLen      =  [this.src.patternVec1(end), this.src.patternVec2(end)]...
                                           * this.src.patternSizeFactor / this.grid.res;
                    switch this.src.plane
                        case 'XY'
                            this.cfg.srcpos    =  [this.src.pos(1:2) -  patternGridLen/2, 0];  
                            this.cfg.srcparam1 = [patternLen(1), 0, 0,  size(this.src.pattern,1)];
                            this.cfg.srcparam2 = [0, patternLen(2), 0,  size(this.src.pattern,2)];
                        case 'YZ'
                            this.cfg.srcpos    = [0, this.src.pos(2:3) - patternGridLen/2];  
                            this.cfg.srcparam1 = [0, patternGridLen(1), 0,  size(this.src.pattern,1)];
                            this.cfg.srcparam2 = [0, 0, patternGridLen(2),  size(this.src.pattern,2)];
                        case 'XZ'
                            this.cfg.srcpos    = [this.src.pos(1)- patternGridLen(1)/2, 0, this.src.pos(3) - patternGridLen(2)/2];  
                            this.cfg.srcparam1 = [patternGridLen(1), 0, 0,  size(this.src.pattern,1)];
                            this.cfg.srcparam2 = [0, 0, patternGridLen(2),  size(this.src.pattern,2)];
                        case 'custom'                            
                            this.cfg.srcparam1=[this.src.planePoint1, size(this.src.pattern,1)]; %srcparam1(1:3)=plane corner, srcparam1(4) pattern Dimension X
                            this.cfg.srcparam2=[this.src.planePoint2, size(this.src.pattern,2)]; %srcparam2(1:3)=plane corner, srcparam2(4) pattern Dimension Y
                    end
                case 'pattern3d'
                    % 'pattern3d' [*] - a 3D illumination pattern. 
                    % srcparam1{x,y,z} defines the dimensions, and srcpattern
                    % is a floating-point pattern array, with values between [0-1].
                    this.cfg.srctype = 'pattern3d';
                case 'fourier'
                    %'fourier' [*] - spatial frequency domain source, similar to 'planar', except
                    % the integer parts of srcparam1(4) and srcparam2(4) represent
                    % the x/y frequencies; the fraction part of srcparam1(4) multiplies
                    % 2*pi represents the phase shift (phi0); 1.0 minus the fraction part of
                    % srcparam2(4) is the modulation depth (M). Put in equations:
                    % S=0.5*[1+M*cos(2*pi*(fx*x+fy*y)+phi0)], (0<=x,y,M<=1)
                    this.cfg.srctype = 'fourier';
                case 'arcsine'
                    % 'arcsine' - similar to isotropic, except the zenith angle is uniform
                    % distribution, rather than a sine distribution.
                    this.cfg.srctype = 'arcsine';
                case 'disk'
                    % 'disk' [*] - a uniform disk source pointing along srcdir; the radius is 
                    % set by srcparam1(1) (in grid unit)
                    this.cfg.srctype = 'disk';
                    this.cfg.srcparam1 = radius;
                case 'fourierx'
                    % 'fourierx' [*] - a general Fourier source, the parameters are 
                    % srcparam1: [v1x,v1y,v1z,|v2|], srcparam2: [kx,ky,phi0,M]
                    % normalized vectors satisfy: srcdir cross v1=v2
                    % the phase shift is phi0*2*pi
                    this.cfg.srctype = 'fourierx';
                case 'fourierx2d'
                    % 'fourierx2d' [*] - a general 2D Fourier basis, parameters
                    % srcparam1: [v1x,v1y,v1z,|v2|], srcparam2: [kx,ky,phix,phiy]
                    % the phase shift is phi{x,y}*2*pi
                    this.cfg.srctype = 'fourierx2d'; 
                case 'zgaussian'
                    %'zgaussian' - an angular gaussian beam, srcparam1(0) 
                    % specifies the variance in the zenith angle
                    this.cfg.srctype = 'zgaussian';
                case 'line'
                    % 'line' - a line source, emitting from the line segment between 
                    % cfg.srcpos and cfg.srcpos+cfg.srcparam(1:3), radiating 
                    % uniformly in the perpendicular direction
                    this.cfg.srctype = 'line';
                    this.cfg.srcparam1 = lineEndPoint;
                case 'slit'
                    % 'slit' [*] - a colimated slit beam emitting from the line segment between 
                    % cfg.srcpos and cfg.srcpos+cfg.srcparam(1:3), with the initial  
                    % dir specified by cfg.srcdir
                    this.cfg.srctype = 'slit';
                    this.cfg.srcparam1 = slitEndPoint;
                case 'pencilarray'
                    % 'pencilarray' - a rectangular array of pencil beams. The srcparam1 and srcparam2
                    % are defined similarly to 'fourier', except that srcparam1(4) and srcparam2(4)
                    % are both integers, denoting the element counts in the x/y dimensions, respectively. 
                    % For exp., srcparam1=[10 0 0 4] and srcparam2[0 20 0 5] represent a 4x5 pencil beam array
                    % spanning 10 grids in the x-axis and 20 grids in the y-axis (5-voxel spacing)
                    this.cfg.srctype = 'pencilarray';
                otherwise
                    error("Unknown source type");
            end
        end
        
        function setOpticalProperties(this)
            %% Parameters
            % cfg.prop: an N by 4 array, each row specifies [mua, mus, g, n] in order.
            %           the first row corresponds to medium type 0
            %           (background) which is typically [0 0 1 1]. The
            %           second row is type 1, and so on. The background
            %           medium (type 0) has special meanings: a photon
            %           terminates when moving from a non-zero to zero voxel.
            % cfg.isreflect: [1]-consider refractive index mismatch, 0-matched index
            % cfg.isnormalized:[1]-normalize the output fluence to unitary source, 0-no reflection
            % cfg.isspecular: 1-calculate specular reflection if source is outside, [0] no specular reflection
            this.cfg.prop(1,:) = [0 0 1 1];
            for i=1:this.op.numOfObjects
                this.cfg.prop(i+1, :) = [this.op.mua(i), this.op.mus(i), this.op.g(i), this.op.ref(i)];
            end
            
            this.cfg.isreflect    = 1;
            this.cfg.isspecular   = 0;
            this.cfg.isnormalized = 1;
        end
        
        function setDetector(this)
            
            
        end
        
        function setSimulation(this)
            %% Parameters:
            % cfg.outputtype: 'flux' - fluence-rate, (default value)
            %           'fluence' - fluence integrated over each time gate, 
            %           'energy' - energy deposit per voxel
            %           'jacobian' or 'wl' - mua Jacobian (replay mode), 
            %           'nscat' or 'wp' - weighted scattering counts for computing Jacobian for mus (replay mode)
            %           for type jacobian/wl/wp
            this.cfg.outputtype   = 'fluence';
            this.cfg.autopilot  = this.mcx.autopilot;
            this.cfg.nblocksize = this.mcx.nblocksize;
            this.cfg.nthread    = this.mcx.nthread;
        end

        function setTime(this)
            %% cfg.tstart:     starting time of the simulation (in seconds)
            % cfg.tstep:      time-gate width of the simulation (in seconds)
            % cfg.tend:       ending time of the simulation (in second) 
            
            this.cfg.tstart = this.time.tstart;
            this.cfg.tstep  = this.time.tstep;
            this.cfg.tend   = this.time.tend;
        end
        
        function phi = calcPhi(this)
             %% Output Parameters:
             %  fluence: a struct array, with a length equals to that of cfg. For each element of fluence:
             %           fluence(i).data: is a 4D array with dimensions specified by [size(vol) total-time-gates]. 
             %                            The content of the array is the normalized fluence at each voxel of each time-gate.
             %           fluence(i).dref: is a 4D array with the same dimension as fluence(i).data
             %                            if cfg.issaveref is set to 1, containing only non-zero values in the 
             %                            layer of voxels immediately next to the non-zero voxels in cfg.vol,
             %                            storing the normalized total diffuse reflectance (summation of the weights 
             %                            of all escaped photon to the background regardless of their direction);
             %                            it is an empty array [] when if cfg.issaveref is 0.
             %           fluence(i).stat: is a structure storing additional information, including
             %                            runtime: total simulation run-time in millisecond
             %                            nphoton: total simulated photon number
             %                            energytot: total initial weight/energy of all launched photons
             %                            energyabs: total absorbed weight/energy of all photons
             %                            normalizer: normalization factor
             %                            unitinmm: same as cfg.unitinmm, voxel edge-length in mm
             %  detphoton: (optional) a struct array, with a length equals to that of cfg.
             %              Starting from v2018, the detphoton contains the below subfields:
             %              detphoton.detid: the ID(>0) of the detector that captures the photon
             %              detphoton.nscat: cummulative scattering event counts in each medium
             %              detphoton.ppath: cummulative path lengths in each medium (partial pathlength)
             %                               one need to multiply cfg.unitinmm with ppath to convert it to mm.
             %              detphoton.mom: cummulative cos_theta for momentum transfer in each medium  
             %              detphoton.p or .v: exit position and direction, when cfg.issaveexit=1
             %              detphoton.w0: photon initial weight at launch time
             %              detphoton.prop: optical properties, a copy of cfg.prop
             %              detphoton.data: a concatenated and transposed array in the order of [detid nscat ppath mom p v w0]'
             %                              "data" is the is the only subfield in all MCXLAB before 2018
             %  vol: (optional) a struct array, each element is a preprocessed volume
             %        corresponding to each instance of cfg. Each volume is a 3D int32 array.
             %  seeds: (optional), if give, mcxlab returns the seeds, in the form of
             %          a byte array (uint8) for each detected photon. The column number
             %          of seed equals that of detphoton.
             %  trajectory: (optional), if given, mcxlab returns the trajectory data for
             %              each simulated photon. The output has 6 rows, the meanings are 
             %              id:  1:    index of the photon packet
             %              pos: 2-4:  x/y/z/ of each trajectory position
             %                5:    current photon packet weight
             %                6:    reserved
             %        By default, mcxlab only records the first 1e7 positions along all
             %        simulated photons; change cfg.maxjumpdebug to define a different limit.
             [fluence, detpt, vol, seeds, traj] = mcxlab(this.cfg);

             if this.time.cw
                 switch this.cfg.outputtype
                     case 'flux'
                        phi = sum(fluence.data,4) *  this.time.tstep;
                     case 'fluence'
                        phi = sum(fluence.data, 4);
                 end
             end
        end
        
    end
end

