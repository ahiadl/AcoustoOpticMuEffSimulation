classdef analysisFunctions
    %ANALYSISFUNCTIONS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods (Static)
        function normMat = normMatf(mat, dim)
            if nargin >1
                numOfDims = length(size(mat));
                dimsOrig = 1:numOfDims;
                dimsNew = dimsOrig;
                dimsNew(dimsNew == dim) = [];
                dimsNew = [dim, dimsNew]; 
                matNew = permute(mat, dimsNew);

                minMat = min(matNew, [], 1);
                maxMat = max(matNew, [], 1);
                span = maxMat - minMat;
                normMat = (matNew - minMat) ./ span;

                dimsRe = [2:dim, 1, dim+1:numOfDims];
                normMat = permute(normMat, dimsRe);
            else
                maxVal = max(mat(:));
                minVal = min(mat(:));
                span = maxVal-minVal;

                normMat = (mat-minVal)/span;
            end
        end 
    end
    
    methods
        function obj = analysisFunctions()
        end
        
    end
end

