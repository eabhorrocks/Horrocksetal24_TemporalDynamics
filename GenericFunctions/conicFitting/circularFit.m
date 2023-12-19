classdef circularFit<ellipticalFit
%A class for executing and managing 2D circle fits
    
    properties (Dependent)
        radius
    end
    
    properties(Dependent,Hidden)
        rad
    end
    
    methods %Accessors
        
        function val=get.radius(obj)
               val=obj.minorDiam/2;
        end

        function obj=set.radius(obj,val)
               obj.majorDiam=2*val;
               obj.minorDiam=2*val;
        end
        
        function val=get.rad(obj)
               val=obj.minorDiam/2;
        end

        function obj=set.rad(obj,val)
               obj.majorDiam=2*val;
               obj.minorDiam=2*val;
        end     
    end
    
    methods
        
        function obj=circularFit(XY)
        %circularFit constructor
        %
        %    obj=circularFit(xy)                    
        %                            
        %IN:                         
        %                            
        %    xy: a 2xN  matrix whose columns are sample coordinates on a circle to be fitted.  
            
            obj@ellipticalFit([]);
            
            if ~nargin||isempty(XY), return; end
            
            obj.XY=XY;
           
            
            XY0=XY;
            [XY,T]=conicFit.homogNorm(XY0);
            
            X=XY(1,:).';
            Y=XY(2,:).';
            
            M= [X.^2 + Y.^2, X, Y, ones(size(X,1),1)];
            
            [~,~,V]=svd(M,0);
            hABC=V(:,end);
            ABC=hABC(2:4).'/hABC(1);
            
            
            if any(~isfinite(ABC))
                warning 'Degenerate circle'
            end
            
            A=ABC(1); B=ABC(2); C=ABC(3);
            
            x0=-A/2; y0=-B/2;
            
            conic=[1 0 -x0 ; 0 1 -y0; -x0 -y0 C];
            conic=T.'*conic*T;
            conic=conic/conic(1);
            
            x0=-conic(end,1);
            y0=-conic(end,2);
            radius=sqrt(x0^2+y0^2-conic(end));
            

            obj.XY=XY0;
            obj.center=[x0,y0];
            obj.majorDiam=2*radius;
            obj.minorDiam=2*radius;
            obj.angle=0;
            
            
        end
        
        function xy = sample(obj,varargin)
        %Generate discrete samples on fitted circle.
        %
        % xy = obj.sample(azimSamps)
        %
        %IN:
        %
        %    azimSamps: azimuthal sample coordinates in degrees.
        %
        %OUT:
        %
        %    xy: A 2x1 cell array of x,y sample locations
        %        on the circle.
            
            assert(numel(varargin)<=1,'Too many input arguments');
            xy = obj.xysim(obj.center,obj.radius,varargin{:});
            xy=num2cell(xy,2);
            
        end       

    end
 
    
    methods (Access = protected, Hidden)
        function propgrp = getPropertyGroups(obj)
            if ~isscalar(obj)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
            else
                propList = struct('center',obj.center,...
                    'radius', obj.rad);
                propgrp = matlab.mixin.util.PropertyGroup(propList);
            end
        end
    end
    
    methods (Static)
        
        
        function obj = groundtruth(xy,center,radius)
        %Creata an object by specifying the geometry instead of fitting it.
        %Useful for plotting a ground truth locus.
        %
        % obj = circularFit.groundtruth(xy,center,radius)
        %
        %IN:
        %
        %        xy: a 2xN matrix of data points (can also be []).
        %    center: ground truth center [x0,y0,z0].
        %    radius:  ground truth circle radius.
     
            obj=circularFit([]);
            obj.XY=xy;
            obj.center=center;
            obj.radius=radius;
            
        end
        
        function xy = xysim(center,radius,varargin)
        % xy = circularFit.xysim(center,R,azimSamps,sig)
        %
        %IN:
        %
        %    center: ground truth circle center [x0,y0]
        %    radius: ground truth circle radius
        %
        %    azimSamps: azimuthal sample coordinates in degrees.
        %    sig: additive Gaussian noise standard deviation.
        %OUT:
        %
        %    xy: A 2xN matrix whose columns are simulated sample locations
        %        on circle.

 
            xy = xysim@ellipticalFit(center,[radius,radius],0,varargin{:});
            
        end
        
    end
end

