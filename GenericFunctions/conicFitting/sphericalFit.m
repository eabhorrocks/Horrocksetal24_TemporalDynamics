classdef sphericalFit<ellipsoidalFit
%A class for executing and managing fits to a sphere
    
    properties (Dependent)
        radius
    end
    
    
    methods %Accessors
        
        function val=get.radius(obj)
               val=obj.aDiam/2;
        end

        function obj=set.radius(obj,val)
               obj.aDiam=2*val;
               obj.bDiam=2*val;
               obj.cDiam=2*val;
        end
        
   
    end
    
    methods
        
        function obj=sphericalFit(XYZ)
        %sphericalFit constructor
        %
        %    obj=sphericalFit(xyz)                    
        %                            
        %IN:                         
        %                            
        %    xyz: a 3xN  matrix whose columns are sample coordinates on a sphere to be fitted.   
        
            obj@ellipsoidalFit([]);
            
            if ~nargin||isempty(XYZ), return; end
            
            obj.XYZ=XYZ;
            
            
            XYZ0=XYZ;
            [XYZ,T]=quadricFit.homogNorm(XYZ0);
            
            
            X=XYZ(1,:).';
            Y=XYZ(2,:).';
            Z=XYZ(3,:).';
            
            M= [X.^2 + Y.^2 + Z.^2, X, Y, Z, ones(size(X,1),1)];
            
            [~,~,V]=svd(M,0);
            hABCD=V(:,end);
            ABCD=hABCD(2:5).'/hABCD(1);
            
            
            if any(~isfinite(ABCD))
                warning 'Degenerate sphere'
            end
            
            A=ABCD(1); B=ABCD(2); C=ABCD(3);D=ABCD(4);
            
            x0=-A/2; y0=-B/2; z0=-C/2;
            
            quadric=diag([1 1 1 D]);
            quadric(1:3,end)=[-x0;-y0;-z0];
            quadric(end,1:3)=[-x0;-y0;-z0].';
            quadric=T.'*quadric*T;
            quadric=quadric/quadric(1);
            
            d=-quadric(1:3,end).';
            R=sqrt(norm(d)^2-quadric(end));
            
            
            obj.center=d;
            obj.aDiam=2*R;
            obj.bDiam=2*R;
            obj.cDiam=2*R;
            obj.ypr=[0,0,0];
            
            
        end
        
        function xyz = sample(obj,varargin)
        %Generate samples on fitted sphere.
        %
        %     xyz = obj.sample(azimSamps,elevSamps)
        %
        %IN:
        %
        %    azimSamps: azimuthal sample coordinates in degrees.
        %    elevSamps: elevation sample coordinates in degrees.
        %
        %OUT:
        %
        %    xyz: A 3x1 cell array of x,y,z sample locations
        %         on surface of the fitted sphere.
            
          assert(numel(varargin)<=2,'Too many input arguments');
          xyz = obj.xyzsim(obj.center,obj.radius,varargin{:});
          xyz=num2cell(xyz,2);
          
        end

    end
 
    
    methods (Access = protected, Hidden)
        function propgrp = getPropertyGroups(obj)
            if ~isscalar(obj)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
            else
                propList = struct('center',obj.center,...
                                  'radius', obj.radius);
                propgrp = matlab.mixin.util.PropertyGroup(propList);
            end
        end
    end
    
    methods (Static)
        
        function obj = groundtruth(xyz,center,radius)
        %Creata an object by specifying the geometry instead of fitting it.
        %Useful for plotting a ground truth surface.
        %
        % obj = sphericalFit.groundtruth(xyz,center,radius)
        %
        %IN:
        %
        %       xyz:  a 3xN matrix of data points (can also be []) 
        %    center:  ground truth center [x0,y0,z0]
        %    radius:  ground truth sphere radius    
        
          
          obj=sphericalFit([]);
          obj.XYZ=xyz;
          obj.center=center;
          obj.radius=radius;
        
        end        
        
        
        
        
        function xyz = xyzsim(center,radius,varargin)
        % xyz = sphericalFit.xysim(center,radius,azimSamps,elevSamps,sig)
        %
        %IN:
        %
        %    center: ground truth center [x0,y0,z0].
        %    radius:  ground truth sphere radius.
        %
        %    azimSamps: azimuthal sample coordinates in degrees.
        %    elevSamps: elevation sample coordinates in degrees.
        %    sig: additive Gaussian noise standard deviation. 
        %
        %OUT:
        %
        %    xyz: A 3xN matrix whose columns are simulated sample locations
        %        on sphere.
            
            xyz =xyzsim@ellipsoidalFit(center,[radius,radius,radius],[0,0,0],varargin{:});
            
        end
        
    end
end

