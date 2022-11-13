classdef illuminationAnalysis<handle
    %ILLUMINATIONANALYSIS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        path
        data
        vars
    end
    
    methods
        function this = illuminationAnalysis()

        end
        
        function [res, vars] = analyse(this, path)
            this.path = path;
            this.load();
            this.processing();
%             this.rotate
%             this.interpolate();
%             this.measure();
            this.plotResults();
            
            res = this.data;
            vars = this.vars;
        end
        
        function load(this)
            tmp = load(this.path);
            this.data.raw = tmp.resCs;
            this.vars = tmp.csVars;
            close(gcf);
            close(gcf);
            
        end
        
        
        function processing(this)
            this.data.maskRaw  = normMatf(mean(this.data.raw, 5));
            mask2 = normMatf(this.data.maskRaw - mean(this.data.maskRaw, 1));
            mask3 = normMatf(mask2 - mean(mask2([1:6, 57:60], :), 1));
            this.data.mask =  mask3.*(mask3>0.1);
            
            this.data.ax1 = this.vars.axDisc1;
            this.data.ax2 = this.vars.axScan;
            
%             this.plotResults();
        end
        
        function plotResults(this)
            figure();
            imagesc(this.vars.disc1Vec, this.vars.scanVec, this.data.mask)
            xlabel(sprintf("%s[mm]", this.vars.axDisc1))
            ylabel(sprintf("%s[mm]", this.vars.axScan))
            colorbar
            axis tight equal
            
%             idx1 = floor(this.vars.scanSizeBin(1)/2);
%             idx2 = floor(this.vars.scanSizeBin(2)/2);
%             
%             figure();
%             subplot(1,2,1)
%             plot(this.vars.disc1Vec, this.data.mask(idx1,:))
%             xlabel(sprintf("%s[mm]", this.vars.axDisc1))
%             ylabel("Illumination [AU]");
%             subplot(1,2,2)
%             plot(this.vars.scanVec, this.data.mask(:,idx2))
%             xlabel(sprintf("%s[mm]", this.vars.axScan))
%             ylabel("Illumination [AU]");
        end
        
        function flush(this)
            this.data = [];
        end
    end
end

