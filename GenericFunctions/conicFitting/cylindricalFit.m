classdef cylindricalFit<quadricFit
%A class for executing and managing cylinder fits
    
    properties (SetAccess=protected)
        center
        radius
        height
    end
    
    properties (Dependent, Hidden)
        
        axis
        yp
        
    end
    
    properties (Hidden)
        
     
        
    end
    
    methods
        
        function ax=get.axis(obj)
            ax=obj.R(:,1);
        end
        
        function yp=get.yp(obj)
            
            yp=[obj.yaw,obj.pitch];
            
        end       
        
        function obj=set.yp(obj,yp)
            
            obj.ypr=[yp,0];
            
        end
        
        function obj=cylindricalFit(XYZ,ax)
        %cylindricalFit class constructor.
        %
        % Syntaxes:
        %
        %    (1) obj=cylindricalFit(xyz)
        %    (2) obj=cylindricalFit(xyz,[y,p])
        %    (3) obj=cylindricalFit(xyz,[ax1,ax2,ax3])
        %                            
        %IN:                         
        %                            
        %        xyz:      A 3xN  matrix whose columns are sample coordinates of a cylinder to be fitted. 
        %      [y,p]:      An initial guess of the yaw and pitch orientation angles (degrees) of the cylinder
        %   [ax1,ax2,ax3]: An initial guess of the cylinder axis 3D direction vector 
        %
        %NOTE: If syntax (1) is used, the routine will auto-generate an
        %      initial guess of the cylinder axis

            if ~nargin, XYZ=[]; end
            
            if isempty(XYZ), return; end
            
            
            XYZ0=XYZ;
            
            [xyzCtr,T]=quadricFit.homogNorm(XYZ);
            %T=eye(4); xyzCtr=XYZ;
            
            if nargin<2
                
                   warnStruct=warning('off');
                   obje=ellipsoidalFit(xyzCtr);
                   warning(warnStruct);
                ax{3}=obje.R(:,1);              
                ax{1}=quadricFit.leastnull(xyzCtr.');
                ax{2}=quadricFit.mostnull(xyzCtr.');

                f(3) = ellipticity(ax{3},xyzCtr);
                f(1) = ellipticity(ax{1},xyzCtr);
                f(2) = ellipticity(ax{2},xyzCtr);

                [~,idx]=min(f);
                ax=ax{idx};
                
                
                params=nan(2,1);
                [params(1),params(2)] = cart2sph(ax(1),ax(2),ax(3));
                
            elseif numel(ax)==2
                
                ax(2)=-ax(2);
                params=ax*pi/180;

            elseif numel(ax)~=3
                
                error 'AX input must be a 2- or 3- element vector ';
                
            end
            
            p0=params;
            
            fun=@(params) ellipticity(params,xyzCtr);
            paramsOpt=fminsearch(fun, p0,optimset('TolFun',1e-10));
            [~, R] = orthproj(xyzCtr,paramsOpt); 
            
            [~,idx]=max(abs(R(:,end)));
            R=sign(R(idx,3))*R;
            R(:,1)=sign(det(R))*R(:,1);
            
            R=R(:,[3,1,2]);
             ypr=quadricFit.rot2ypr(R);
             ypr(3)=0;
            R=quadricFit.ypr2rot( ypr );
            
            s1=sign(R(1)); if ~s1, s1=+1; end
            s2=sign(R(5)); if ~s2, s2=+1; end
            R=R.*[s1,s2,1];
            R(:,3)=cross(R(:,1),R(:,2));
            
            xyzCtr=R.'*xyzCtr;
 
            objCirc=circularFit(xyzCtr(2:3,:));
            
            xmin=min(xyzCtr(1,:));
            xmax=max(xyzCtr(1,:));
            
            hcenter=[0; objCirc.center(:);1];
            
            H=T\blkdiag(R,1);
            center=H*hcenter; center=center(1:3)/center(4);
            
            
            obj.XYZ=XYZ0;
            obj.center=center(:).';
            obj.radius=det(H)^(1/3)*objCirc.radius;
            obj.ypr=quadricFit.rot2ypr(R);
            obj.height=det(H)^(1/3)*(xmax-xmin);
            
        end
        
        function xyz = sample(obj,varargin)
        %Generate discrete samples on fitted cylinder.
        %
        % xyz = obj.sample(azimSamps,elevSamps)
        %
        %IN:
        %
        %    azimSamps:  Azimuthal sample coordinates in degrees.
        %    axialSamps: Axial coordinates of samples (along axis relative to center).
        %                These values should satisfy max(axialSamps)=-min(axialSamps),
        %                otherwise the code will re-center them. 
        %
        %OUT:
        %
        %    xyz: A 3x1 cell array of x,y,z sample locations
        %         on the fitted cylinder.
            
          assert(numel(varargin)<=2,'Too many input arguments');
          xyz = obj.xyzsim(obj.center,obj.radius,obj.yp,varargin{:});
          xyz=num2cell(xyz,2);  
          
        end
        
        function varargout=fimplicit(obj,varargin)
            
            varargin=quadricFit.surfDefaults(varargin);
            
            x0=obj.center(1);
            y0=obj.center(2);
            z0=obj.center(3);
            rad=obj.radius;
            Rt=obj.R.';
            hhalf=obj.height/2;
            
            Q=diag([0,1,1]./rad^2);
            
            hFit=fimplicit3(@quadform,varargin{:});
            
            if nargout, varargout={hFit}; end
                        xlabel 'X'; ylabel 'Y', zlabel 'Z'
            function d=quadform(x,y,z)
                
                xyz=Rt*[x(:).'-x0; y(:).'-y0; z(:).'-z0];
                
                outrange =  abs( xyz(1,:) )>hhalf ;

                d=reshape( (sum((Q*xyz).*xyz)-1), size(x));
                
                d(outrange)=nan;
            end
            
            
        end
        
        
        function varargout=surf(obj,varargin)
            
            varargin=quadricFit.surfDefaults(varargin);
            
            [X,Y,Z]=cylinder(obj.radius,1e3); 
            Z=(Z-0.5)*obj.height;
            
            T0=[ 0     0     1
                 0     1     0
                -1     0     0];
            
            T=[obj.R*T0,obj.center(:); 0 0 0 1];
            
            
            hFit=surf(X,Y,Z, 'Parent',hgtransform('Matrix',T),varargin{:});
            
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
                                  'radius', obj.radius,...
                                  'height', obj.height,...
                                  'yaw',obj.yaw,...
                                  'pitch',obj.pitch);
                propgrp = matlab.mixin.util.PropertyGroup(propList);
            end
        end
        

        
    end
    
    methods (Static)
        
        function obj = groundtruth(xyz,center,radius,height,yp)
        %Creata an object by specifying the geometry instead of fitting it.
        %Useful for plotting a ground truth surface.
        %
        % obj = cylindricalFit.groundtruth(xyz,center,radius,yp)
        %
        %IN:
        %
        %       xyz:  a 3xN matrix of data points (can also be []) 
        %    center:  ground truth center [x0,y0,z0]
        %    radius:  ground truth cylinder radius
        %    height:  ground truth cylinder height
        %        yp:  ground truth [yaw,pitch] of the cylinder in degrees
        
          
          obj=cylindricalFit([]);
          obj.XYZ=xyz;
          obj.center=center;
          obj.radius=radius;
          obj.height=height;
          obj.yp=yp;
                  
        end        
        
        
        
        function xyz = xyzsim(center,radius,yp,theta,h,sig)
        % xyz = cylindricalFit.xysim(center,radii,yp,azimSamps,elevSamps,sig)
        %
        %IN:
        %
        %    center:  ground truth cylinder center [x0,y0,z0]
        %    radius:  ground truth radius
        %        yp:  yaw, pitch angle vector [y,p] in degrees.
        %
        %    azimSamps:  Azimuthal sample coordinates in degrees.
        %    axialSamps: Axial coordinates of samples (along axis relative to center).
        %                These values should satisfy max(axialSamps)=-min(axialSamps),
        %                otherwise the code will re-center them. 
        %    sig: additive Gaussian noise standard deviation. 
        %
        %OUT:
        %
        %    xyz: A 3xN matrix whose columns are simulated sample locations
        %         on cylinder.
            
            
            if ~exist('sig','var')||isempty(sig), sig=0; end
            if ~exist('theta','var')||isempty(theta), theta=0:10:350; end
            if ~exist('h','var')||isempty(h), h=linspace(-2*radius,2*radius,20); end
            
            [theta,h]=ndgrid(theta,h-(max(h)+min(h))/2);
            
            [X,Y,Z]=pol2cart(theta*pi/180,radius,h);
            

            xyz=[X(:),Y(:),Z(:)].';
            
            
            if sig
                
                N=size(xyz,2);
                %dr=normalize(xyz.*[1;1;0],1,'norm').*randn(1,N)*sig;
                dr=normalize(randn(3,N),1,'norm')*sig;
                
                xyz=xyz+dr; %add radial noise
                
            end
            
            xyz=xyz([3,1,2],:);
            
            
            R=quadricFit.ypr2rot([yp,0]);
            
            xyz=R*xyz+center(:);
            
        end
        
    end
end

%%%CLASS-RELATED FCNs

function E = ellipticity(params,XYZ)

  xy=orthproj(XYZ,params);
  
  try
      obj=ellipticalFit(xy);
      dists=vecnorm(xy-obj.center(:));
      E=var(dists);
  catch
     E=inf;
  end
  
  
end

function  [xy,basis] = orthproj(XYZ,params)


    n=numel(params);
    if n==2

        ax=nan(3,1);
        [ax(1),ax(2),ax(3)] = sph2cart(params(1),params(2),1);

    elseif n==3

        ax=params(:);

    end

    basis=null(ax(:).');

    xy=basis.'*XYZ;

    if nargout==2
        basis=[basis,ax];
    end
  
  
end

