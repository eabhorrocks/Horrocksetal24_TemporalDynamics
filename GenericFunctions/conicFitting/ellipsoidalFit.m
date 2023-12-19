classdef ellipsoidalFit<quadricFit
%A class for executing and managing ellipsoid fits
    
    properties (SetAccess=protected, Hidden)
        
        aDiam
        bDiam 
        cDiam
        

        
    end
    
    properties (SetAccess=protected)
        center 
    end
    
    properties (Dependent,Hidden) 
        abc 
    end
 
    properties (Dependent) %Not settable
        a  
        b 
        c
    end

    
    methods %Accessors
        
        function val=get.a(obj)
               val=obj.aDiam/2;
        end   
        function val=get.b(obj)
               val=obj.bDiam/2;
        end
        function val=get.c(obj)
               val=obj.cDiam/2;
        end
     
        function obj=set.abc(obj,val)
            
               val=2*sort(val,'descend');
               obj.aDiam=val(1);
               obj.bDiam=val(2);
               obj.cDiam=val(3);              
        end
        
        function val=get.abc(obj)
            
               val=[obj.aDiam, obj.bDiam, obj.cDiam]/2;
        end

        
        
    end
    
    methods
        
        function obj=ellipsoidalFit(XYZ)
        %ellipsoidalFit constructor
        %
        %    obj=ellipseFit(xyz)                    
        %                            
        %IN:                         
        %                            
        %    xyz: a 3xN  matrix whose columns are sample coordinates on an ellipsoid to be fitted.   
            
            
            if ~nargin, XYZ=[]; end

            if isempty(XYZ), return; end
            
            
            XYZ0=XYZ;
            [XYZ,T]=quadricFit.homogNorm(XYZ0);
            
            X=XYZ(1,:).';
            Y=XYZ(2,:).';
            Z=XYZ(3,:).';
            
            e=+ones(size(X,1),1);
            M= [X.^2, X.*Y, X.*Z, X, ...
                      Y.^2, Y.*Z, Y, ...
                            Z.^2  Z, ...
                                  e];
            
           
            ABCDEFGHIJ=quadricFit.mostnull(M);
            ABCDEFGHIJ=num2cell(ABCDEFGHIJ);
 
            [A,B,C,D,E,F,G,H,I,J]=deal(ABCDEFGHIJ{:});
            
            
            Q=[A, B, C; %D
               0  E, F; %G
               0  0  H];%I
                        %J

            Q=Q/2+Q.'/2;            
                        
            W=T.'*[Q,[D;G;I]/2;[D,G,I]/2,J]*T;
            
            Q=W(1:3,1:3);
            x0=-Q\W(1:3,end);
            
            T=eye(4); T(1:3,4)=x0;
            
            W=T.'*W*T; W=-W/W(end);
            
            [R,eigen]=eig(W(1:3,1:3),'vector');
            
            if ~all(eigen>=0), warning 'Data is too noisy or non-ellipsoidal. Poor quality fit is anticipated.'; end
            idx=eigen<=0;
            eigen=abs(eigen);
            eigen(idx)=max(eigen)/1000;

            AxesDiams = 2*sqrt(1./eigen);
            
            [AxesDiams,idx] = sort(AxesDiams(:).','descend');
            R=R(:,idx);
            
            s1=sign(R(1)); if ~s1, s1=+1; end
            s2=sign(R(5)); if ~s2, s2=+1; end
            R=R.*[s1,s2,1];
            R(:,3)=cross(R(:,1),R(:,2));
            
            
            obj.XYZ=XYZ0;
            obj.center=x0(:).';
            obj.aDiam=AxesDiams(1);
            obj.bDiam=AxesDiams(2);
            obj.cDiam=AxesDiams(3);
            obj.R=R;
            
            
        end
        
        function xyz = sample(obj,varargin)
        %Generate samples on fitted ellipsoid.
        %
        % xyz = obj.sample(azimSamps,elevSamps)
        %
        %IN:
        %
        %    azimSamps: azimuthal sample coordinates in degrees.
        %    elevSamps: elevation sample coordinates in degrees.
        %
        %OUT:
        %
        %    xyz: A 3x1 cell array of x,y,z sample locations
        %        on the fitted ellipsoid surface.
            
          assert(numel(varargin)<=2,'Too many input arguments');
          xyz = obj.xyzsim(obj.center,obj.abc,obj.ypr,varargin{:});
          xyz=num2cell(xyz,2);
            
        end
        

        function varargout=fimplicit(obj,varargin)
            
            varargin=quadricFit.surfDefaults(varargin); 
            
            x0=obj.center(1);
            y0=obj.center(2);
            z0=obj.center(3);
            a=obj.a;
            b=obj.b;
            c=obj.c;
            R=obj.R;

            Q=R*diag([1./a^2;1./b^2;1./c^2])*R.';

            hFit=fimplicit3(@quadform,varargin{:});
           
            if nargout, varargout={hFit}; end
                        
            function d=quadform(x,y,z)

                xyz=[x(:).'-x0; y(:).'-y0; z(:).'-z0];
                
                d=reshape( sum((Q*xyz).*xyz)-1 , size(x));
                
            end
            
            
        end
        
        function varargout=surf(obj,varargin)
            
            varargin=quadricFit.surfDefaults(varargin);
            
            [X,Y,Z]=sphere(1e2);
            
            T=[obj.abc.*obj.R,obj.center(:); 0 0 0 1];
            
            
            hFit=surf(X,Y,Z, 'Parent',hgtransform('Matrix',T),varargin{:});
                        xlabel 'X'; ylabel 'Y', zlabel 'Z'
            if nargout, varargout={hFit}; end
            
        end

        function varargout= showfit(obj,varargin)
            
            [varargout{1:nargout}]=obj.surf(varargin{:});
            showfit@quadricFit(obj);
        end      
        
    end
    
    
    methods (Access = protected, Hidden)
        function propgrp = getPropertyGroups(obj)
            if ~isscalar(obj)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
            else
                propList = struct('center',obj.center,...
                    'a',obj.a,...
                    'b',obj.b,...
                    'c',obj.c,...
                    'yaw',obj.yaw,...
                    'pitch',obj.pitch,...
                    'roll',obj.roll);
                propgrp = matlab.mixin.util.PropertyGroup(propList);
            end
        end
        

    end

    
    methods (Static)
        
        function obj = groundtruth(xyz,center,abc,ypr)
        %Creata an object by specifying the geometry instead of fitting it.
        %Useful for plotting a ground truth surface.
        %
        % obj = ellipsoidalFit.groundtruth(xyz,center,abc,ypr)
        %
        %IN:
        %
        %      xyz: a 3xN matrix of data points (can also be []) 
        %    center: ground truth center [x0,y0,z0]
        %    radii:  ground truth axis radii [a,b,c] in descending order
        %    ypr: yaw, pitch, roll angle vector [y,p,r] in degrees.      
        
          
          obj=ellipsoidalFit([]);
          obj.XYZ=xyz;
          obj.center=center;
          obj.abc=abc;
          obj.ypr=ypr;
        
        end
        
        function xyz = xyzsim(center,abc,ypr,azim,elev,sig)
        % xyz = ellipsoidalFit.xysim(center,radii,ypr,azimSamps,elevSamps,sig)
        %
        %IN:
        %
        %    center: ground truth center [x0,y0,z0]
        %    radii:  ground truth axis radii [a,b,c] in descending order
        %    ypr: yaw, pitch, roll angle vector [y,p,r] in degrees.
        %
        %    azimSamps: azimuthal sample coordinates in degrees.
        %    elevSamps: elevation sample coordinates in degrees.
        %    sig: additive Gaussian noise standard deviation. 
        %
        %OUT:
        %
        %    xyz: A 3xN matrix whose columns are simulated sample locations
        %        on ellipsoid.
            
            assert(issorted(abc,'descend'),"Ellipse diameters abc must be in descending order")
            
            if ~exist('sig','var')||isempty(sig), sig=0; end
            if ~exist('azim','var')||isempty(azim), azim=0:10:350; end
            if ~exist('elev','var')||isempty(elev), elev=0:10:350; end
            
            [azim,elev]=ndgrid(azim,elev);
            
            [X,Y,Z]=sph2cart(azim*pi/180,elev*pi/180,1);
            

            xyz=([X(:),Y(:),Z(:)].*abc).';
            
            
            if sig
                
                N=size(xyz,2);
               %dr=normalize(xyz,1,'norm').*randn(1,N)*sig;
                dr=normalize(randn(3,N),1,'norm')*sig;
                
                xyz=xyz+dr; %add radial noise
                
            end
            
            R=quadricFit.ypr2rot(ypr);
            
            xyz=R*xyz+center(:);
            
        end
        
    end
    
end

