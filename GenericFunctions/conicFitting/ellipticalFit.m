classdef ellipticalFit<conicFit
%A class for executing and managing 2D ellipse fits
    
    properties (SetAccess=protected, Hidden)
        
        majorDiam
        minorDiam 
 
    end
    
    properties (SetAccess=protected)
        center 
    end
    
    properties (Dependent,Hidden) 
        ab 
    end
 
    properties (Dependent) %Not settable
        a  
        b  
    end

    
    methods %Accessors
        
        function val=get.a(obj)
               val=obj.majorDiam/2;
        end
        
        function val=get.b(obj)
               val=obj.minorDiam/2;
        end


     
        function obj=set.ab(obj,val)
            
               val=2*sort(val,'descend');
               obj.majorDiam=val(1);
               obj.minorDiam=val(2);
        end
        
        function val=get.ab(obj)
            
               val=[obj.majorDiam, obj.minorDiam]/2;
        end

        
        
    end
    
    methods
        
        function obj=ellipticalFit(XY)
        %ellipticalFit constructor
        %
        %    obj=ellipseFit(xy)                    
        %                            
        %IN:                         
        %                            
        %    xy: a 2xN  matrix whose columns are sample coordinates on an ellipse to be fitted.    
            
            
            if ~nargin, XY=[]; end

            if isempty(XY), return; end
            
            
            XY0=XY;
            [XY,T]=conicFit.homogNorm(XY0);
            
            X=XY(1,:).';
            Y=XY(2,:).';
            
            M= [X.^2, X.*Y, Y.^2, X, Y, +ones(size(X,1),1)];
            
            
            
            ABCDEF=conicFit.mostnull(M);
            ABCDEF=num2cell(ABCDEF);          
            [A,B,C,D,E,F]=deal(ABCDEF{:});
            
            
            Q=[A, B/2;B/2 C];
            
            J=T.'*[Q,[D;E]/2;[D,E]/2,F]*T;
            
            [A,B,D,C,E,F]=deal(J(1),2*J(2),2*J(3),J(5),2*J(6),-J(9));
            
            Q=[A, B/2;B/2 C];
            x0=-Q\[D;E]/2;
            
            
            dd=F+x0'*Q*x0;
            
            Q=Q/dd;
            
            [R,eigen]=eig(Q);
            eigen=eigen([1,4]);
            
            if ~all(eigen>=0), warning 'Data is too noisy or non-elliptical. Poor quality fit is anticipated.'; end
            idx=eigen<=0;
            eigen=abs(eigen);
            eigen(idx)=max(eigen)/1000;

            AxesDiams = 2*sqrt(1./eigen);
            
            [AxesDiams,idx] = sort(AxesDiams(:).','descend');
            R=R(:,idx);
            
            e=R(:,1);
            
            if ~e(1)
                theta=90;
            else
                e=e*sign(e(1));
                theta=atan2d(e(2)+eps,e(1));
            end
            
            obj.XY=XY0;
            obj.center=x0(:).';
            obj.majorDiam=AxesDiams(1);
            obj.minorDiam=AxesDiams(2);
            obj.angle=theta;
            
            
        end
        
        function xy = sample(obj,varargin)
        %Generate discrete samples on the fitted ellipse.
        %
        % xy = sample(obj,azimSamps)
        %
        %IN:
        %
        %    azimSamps: azimuthal sample coordinates in degrees.
        %
        %OUT:
        %
        %    xy: A 2x1 cell array of sample locations
        %        on the fitted ellipse.
            
          assert(numel(varargin)<=1,'Too many input arguments');
          xy = obj.xysim(obj.center,obj.ab,obj.theta,varargin{:});
          xy=num2cell(xy,2);
             
        end
        

        function hFit=fimplicit(obj,varargin)
            
            varargin=conicFit.lineDefaults(varargin);
            
            x0=obj.center(1);
            y0=obj.center(2);
            a=obj.a;
            b=obj.b;
            theta=obj.theta;
            
            R2d = @(x) [cosd(x),-sind(x); sind(x),cosd(x)];
            R=R2d(-theta);
            Q=R.'*diag([1./a^2;1./b^2])*R;

            hFit=fimplicit(@quadform,varargin{:});
           
            
            function z=quadform(x,y)

                xy=[x(:).'-x0; y(:).'-y0];
                
                z=reshape( sum((Q*xy).*xy)-1 , size(x));
                
            end
            
            
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
                    'angle',obj.theta);
                propgrp = matlab.mixin.util.PropertyGroup(propList);
            end
        end
    end

    
    methods (Static)
        
        function obj = groundtruth(xy,center,ab, theta)
        %Creata an object by specifying the geometry instead of fitting it.
        %Useful for plotting a ground truth locus.
        %
        % obj = ellipticFit.groundtruth(xy,center,ab, theta)
        %
        %IN:
        %
        %        xy: a 2xN matrix of data points (can also be [])
        %    center: ground truth center [x0,y0,z0]
        %    ab:  ground truth axis radii [a,b] in descending order
        %    theta: ellipse rotation angle in degrees.
            
            
            obj=ellipticalFit([]);
            obj.XY=xy;
            obj.center=center;
            obj.ab=ab;
            obj.angle=theta;
            
        end
        
        function xy = xysim(center,ab,theta,t,sig)
        % xy = ellipticalFit.xysim(center,radii,angle,azimSamps,sig)
        %
        %IN:
        %
        %    center: ground truth ellipse center [x0,y0]
        %    radii: ground truth major/minor radii [a,b]
        %    angle: ground truth ellipse tilt angle in degrees
        %
        %    azimSamps: azimuthal sample coordinates in degrees.
        %    sig: additive Gaussian noise standard deviation. 
        %
        %OUT:
        %
        %    xy: A 2xN matrix whose columns are simulated sample locations
        %        on ellipse.
        
            center=num2cell(center); ab=num2cell(ab);
            [x0,y0]=deal(center{:});
            [a,b]=deal(ab{:});
           
        
            assert(a>=b,"a must be greater than or equal to b")
            
            if ~exist('sig','var')||isempty(sig), sig=0; end
            if ~exist('t','var')||isempty(t), t=0:10:350; end
            
            R2d = @(x) [cosd(x),-sind(x); sind(x),cosd(x)];
            fun=@(t) bsxfun(@plus,R2d(theta)*[a*cosd(t+theta);b*sind(t+theta)],[x0;y0]);
            
            xy=fun(t);
            
            if sig
                
                N=numel(t);
                %dr=normalize(xy-[x0;y0],1,'norm').*randn(1,N)*sig;
                dr=normalize(randn(2,N),1,'norm')*sig;
                
                xy=xy+dr; %add radial noise
                
            end
            
            
        end
        
    end
    
end

