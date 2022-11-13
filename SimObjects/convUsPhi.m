classdef convUsPhi < handle
    %CONVUSPHI Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vars
        us
        phi
    end
    
    methods (Static)
        function uVars = createUserVars() 
           uVars.usAx      = [];
           uVars.usDir     = [];
           uVars.usFocalDist = [];
        end
    end
    
    methods
        function this = convUsPhi()
            
        end
        
        function setVars(this, uVars)
            this.vars.usAx = uVars.usAx;
            this.vars.usDir = uVars.usDir;
            this.vars.usFocalDist = uVars.usFocalDist; % from Boundary;
            
            this.vars.xVecGrid = uVars.xVecGrid;
            this.vars.yVecGrid = uVars.yVecGrid;
            this.vars.zVecGrid = uVars.zVecGrid;
            
            this.vars.xLen = length(this.vars.xVecGrid);
            this.vars.yLen = length(this.vars.yVecGrid);
            this.vars.zLen = length(this.vars.zVecGrid);
            
            switch this.vars.usAx
                case 'X'
                    this.vars.axVecRaw  = this.vars.xVecGrid; 
                    this.vars.tr1VecRaw = this.vars.yVecGrid;
                    this.vars.tr2VecRaw = this.vars.zVecGrid;
                    
                    this.vars.dimAlignVec = [2,3,1];

                case 'Y'
                    this.vars.axVecRaw  = this.vars.yVecGrid; 
                    this.vars.tr1VecRaw = this.vars.xVecGrid;
                    this.vars.tr2VecRaw = this.vars.zVecGrid;

                    this.vars.dimAlignVec = [1,3,2];
                    
                case 'Z'
                    this.vars.axVecRaw  = this.vars.zVecGrid; 
                    this.vars.tr1VecRaw = this.vars.yVecGrid;
                    this.vars.tr2VecRaw = this.vars.xVecGrid;

                    this.vars.dimAlignVec = [2,1,3];
            end
            
            % Matching resolution of Phi(low) to US(high):
            this.vars.axVec  = this.vars.axVecRaw(1)  : this.vars.us.dAx  : this.vars.axVecRaw(end);
            this.vars.tr1Vec = this.vars.tr1VecRaw(1) : this.vars.us.dtr1 : this.vars.tr1VecRaw(end);
            this.vars.tr2Vec = this.vars.tr2VecRaw(1) : this.vars.us.dtr2 : this.vars.tr2VecRaw(end);
            
            this.vars.axSize  = length(this.vars.axVec);
            this.vars.tr1Size = length(this.vars.tr1Vec);
            this.vars.tr2Size = length(this.vars.tr2Vec);
            
            this.alignUS();
            this.calcPadding();
            this.calcConvIdxs();
        end
        
        function vars = getVars(this)
           vars = this.vars;
           vars.us  = [];
           vars.phi = [];
        end
        
        function setUSData(this, usData, usVars)
            this.us.pulses = usData.pulses;
            this.vars.us = usVars;
            
            this.vars.tr1USLenIdx = size(usData.pulses, 1);
            this.vars.tr2USLenIdx = size(usData.pulses, 2);
            this.vars.pulseSizeIdx = this.vars.us.pulseSizeIdx;
        end
        
        function alignUS(this)
            fprintf("Aligning US Pulses to Phi Grid\n");
            if this.vars.usDir == 1
                this.us.pulsesAligned = this.us.pulses;
                this.vars.pulsDistFromTransAligned = this.vars.us.depthAxPulses;
                this.vars.pulseDistAxAligned       = this.vars.pulsDistFromTransAligned - this.vars.us.focalLen + this.vars.usFocalDist;
            elseif this.vars.usDir == -1
                this.us.pulsesAligned = flip(this.us.pulses,4);
                this.vars.pulsDistFromTransAligned = flip(this.vars.us.depthAxPulses);
                this.vars.pulseDistAxAligned       = -(this.vars.pulsDistFromTransAligned - this.vars.us.focalLen) + this.vars.usFocalDist;
            end
            [gridsOffset, i] = min(abs(this.vars.pulseDistAxAligned));
            this.vars.pulseDistAxAligned = this.vars.pulseDistAxAligned - gridsOffset * sign(this.vars.pulseDistAxAligned(i));
        end
        
        function calcPadding(this)
            fprintf("Calculate Padded Matrix Size for Convolution\n");
            % Calculate Padded Matrix Size:
            this.vars.tr1PadSize = this.vars.tr1USLenIdx-1;
            this.vars.tr2PadSize = this.vars.tr2USLenIdx-1;
            this.vars.axPadSize = 2*(this.vars.pulseSizeIdx-1);
            
            this.vars.phiBB = [ this.vars.tr1PadSize/2+1,   this.vars.tr2PadSize/2+1, this.vars.pulseSizeIdx;...
                                this.vars.tr1PadSize/2+this.vars.tr1Size,   this.vars.tr2PadSize/2+this.vars.tr2Size, this.vars.pulseSizeIdx-1 + this.vars.axSize];
            
            this.vars.padedSize = [this.vars.tr1Size,    this.vars.tr2Size,    this.vars.axSize] + ...
                                  [this.vars.tr1PadSize, this.vars.tr2PadSize, this.vars.axPadSize];
            
            this.vars.axVecConv = ((1:1:this.vars.phiBB(2,3))- this.vars.phiBB(1,3))* this.vars.us.dAx;
        end
        
        function calcConvIdxs(this)
            fprintf("Calculating convolution indices\n");
            this.vars.overlapIdx = find(this.vars.pulseDistAxAligned == 0);
            this.vars.endIdx     = this.vars.overlapIdx + this.vars.axSize - 1;
            this.vars.startIdx   = this.vars.overlapIdx - (this.vars.pulseSizeIdx-2);
        end
        
        function setPhi(this, phiData, phiVars)
            fprintf("Interpolating phi to match US\n");
            this.phi.raw = phiData;
            this.vars.phi = phiVars;
            
            this.phi.aligned = permute(this.phi.raw, this.vars.dimAlignVec);
            
            [X, Y, Z]    = meshgrid(this.vars.tr1VecRaw, this.vars.tr2VecRaw, this.vars.axVecRaw);
            [Xq, Yq, Zq] = meshgrid(this.vars.tr1Vec, this.vars.tr2Vec, this.vars.axVec);
            this.phi.int = normMatf(interp3(X, Y, Z, this.phi.aligned, Xq, Yq, Zq, 'cubic'));

            this.phi.pad = zeros(this.vars.padedSize);
            
            this.phi.pad(this.vars.phiBB(1,1):this.vars.phiBB(2,1), ...
                         this.vars.phiBB(1,2):this.vars.phiBB(2,2), ...
                         this.vars.phiBB(1,3):this.vars.phiBB(2,3)) = this.phi.int;

        end

        function conv3DGPU(this)
           fprintf("Convolving\n");
           reset(gpuDevice);
           this.phi.conv = zeros(size(this.phi.int,1), size(this.phi.int,2), this.vars.phiBB(2,3));
           phiPad  = gpuArray(this.phi.pad);
           tic
           for k = 1:this.vars.phiBB(2,3)
               if ~mod(k,100); fprintf("%d\n", k); end
               curPulseIdx =  this.vars.startIdx+k-1;
               pulse    = gpuArray(squeeze(this.us.pulsesAligned(:,:,curPulseIdx,:)));
               pulseMat = gpuArray(zeros(this.vars.tr1USLenIdx, this.vars.tr1USLenIdx, this.vars.padedSize(3)));
               pulseMat(:, :, k:k+size(pulse,3)-1) = pulse;
               this.phi.conv(:,:,k) = gather(convn(phiPad, flip(flip(flip(pulseMat,1),2),3),'valid'));
           end
           toc
           fprintf("Done Convolving!\n");
        end

        function res = calcConv(this, phi, phiVars)
            this.setPhi(phi, phiVars);
            this.conv3DGPU();
            res = this.phi.conv;
        end
        
    end
end

%         function conv(this)
%            fprintf("Convolving\n");
%            convMat = zeros(size(this.phi.int,1), size(this.phi.int,2), this.vars.padedSize(3));
%            phiPad = this.phi.pad;
%            tic
%            for k = 1:this.vars.padedSize(3)
%                curPulseIdx =  this.vars.startIdx+k-1;
%                pulse = squeeze(this.us.pulsesAligned(:,:,curPulseIdx,:));
%                pulseMat = zeros(this.vars.padedSize);
%                pulseMat(1:this.vars.tr1USLenIdx, 1:this.vars.tr1USLenIdx, k:k+size(pulse,3)-1) = pulse;
%                for i = 1:size(this.phi.int,1)
%                    pulseMatConv = circshift(pulseMat, i, 1);
%                    if (i==1 && k==1) ;tic; end
%                    parfor j = 1:size(this.phi.int,2)
%                         pulseMatConv2 = circshift(pulseMatConv, j, 2);
%                         res(j) = sum(sum(sum(pulseMatConv2.*phiPad,1),2),3);  
%                    end
%                    if (i==1 && k==1) ;toc; end
%                    convMat(i,:,k) = res;
%                end
%            end
%            toc
%            fprintf("Done Convolving!\n");
%         end

%            figure(); 
%            subplot(2,2,1)
%            imagesc(squeeze(pulseMat(:,:,1)))
%            colorbar
%            subplot(2,2,2)
%            imagesc(squeeze(pulseMat(:,:,116)))
%            colorbar
%            subplot(2,2,3)
%            imagesc(squeeze(this.phi.conv(:,129,:)))
%            colorbar
%            subplot(2,2,4)
%            imagesc(squeeze(this.phi.conv(:,:,1)))
%            colorbar
%            
%            figure();
%            plot(squeeze(pulseMat(9,9,:))); hold on
%            yyaxis right
%            plot(squeeze(phiPad(129,129,:)))
%            
%            close all
%            figure(); 
%            subplot(2,2,1)
%            imagesc(squeeze(pulseMat(:,:,end)))
%            colorbar
%            subplot(2,2,2)
%            imagesc(squeeze(this.phi.conv(:,:,1)))
%            xlim([113, 129])
%            ylim([113, 129])
%            colorbar
%            subplot(2,2,3)
%            plot(squeeze(pulseMat(:,9,:))')
%            subplot(2,2,4)
%            imagesc(squeeze(phiPad(:,129,:)))