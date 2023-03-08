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
            mask3 = normMatf(mask2 - mean(mask2(1:6, :), 1));
            this.data.mask =  mask3.*(mask3>0.1);
            
            this.data.ax1 = this.vars.axDisc1;
            this.data.ax2 = this.vars.axScan;

        end
        
        function plotResults(this)
            ax1 = this.vars.disc1Vec - mean(this.vars.disc1Vec);
            ax2 = this.vars.scanVec - mean(this.vars.scanVec);

            figure();
            subplot(1,2,1)
            imagesc(ax1, ax2, this.data.mask)
            xlabel(sprintf("%s[mm]", this.vars.axDisc1))
            ylabel(sprintf("%s[mm]", this.vars.axScan))
            h = colorbar;
            ylabel(h, "Linear[AU]")
            axis tight equal
            subplot(1,2,2)
            imagesc(ax1, ax2, db(this.data.mask))
            xlabel(sprintf("%s[mm]", this.vars.axDisc1))
            ylabel(sprintf("%s[mm]", this.vars.axScan))
            h = colorbar;
            ylabel(h, "dB");
            axis tight equal

        end
        
        function flush(this)
            this.data = [];
        end
    end
end

