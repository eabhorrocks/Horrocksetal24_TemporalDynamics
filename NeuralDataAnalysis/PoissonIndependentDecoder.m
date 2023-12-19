classdef PoissonIndependentDecoder
    
    % Object to perform basic decoding
    
    properties
        tuningAxes; % cell array of tuning axes
        tuningDims;
        tuning; % struct array (cells) with tuning field
        allTuning
        bias; % sum of all tuning for all cells
        
    end
    
    methods
        function obj = PoissonIndependentDecoder(varargin) %
        end
        
        function obj = trainDecoder(obj, tuningStructArray)
            % assign tuning array struct to properties
            obj.tuning = tuningStructArray;
            
            % calculate the bias (sum of tuning curves)
            obj.tuningDims = numel(size(tuningStructArray(1).tuning));
            obj.allTuning = cat(obj.tuningDims+1,tuningStructArray.tuning);
            
            obj.bias = sum(obj.allTuning,obj.tuningDims+1);
        end
        
        function pred = testDecoder(obj, spikeCounts)

            if numel(spikeCounts)~=numel(obj.tuning)
                error('must provide spike counts for each cell trained on')
            end
            
            % initialise loglikelihood distribution
            loglikelihood = zeros(size(obj.tuning(1).tuning));

            for iunit = 1:numel(spikeCounts) % loop through cells in population
                loglikelihood = loglikelihood +...
                    spikeCounts(iunit)*log(obj.tuning(iunit).tuning);
            end

            % slower than for loop above
            %logg2 = sum(bsxfun(@times, log(obj.allTuning), reshape(spikeCounts,1,1,length(spikeCounts))),3);
            
            % account for bias due to inhomegenous tuning
            loglikelihood = loglikelihood - obj.bias;
            
           % imagesc(loglikelihood)
           % colormap(bone), colorbar
            [~, idx] = max(loglikelihood(:)); % find columnwise index of max
            % convert to ndim co-ordinates
            numdims = sum(size(loglikelihood)>1);
            [coords(1:numdims).idx] = ind2sub(size(loglikelihood),idx);
            
            for idim = 1:numdims
                pred(idim) = obj.tuningAxes{idim}(coords(idim).idx);
            end
            
            
        end
        
    end
    
end
