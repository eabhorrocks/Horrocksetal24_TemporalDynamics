classdef rightcircularconeFit<quadricFit
%A class for executing and managing right circular cone fits
    
    properties (SetAccess=protected)
        vertex (1,3)
        cone_angle (1,1)
        height (1,1) %upper and lower extremes along axis of the cone
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
        
        
        
        function obj=rightcircularconeFit(XYZ,varargin)
        %rightcircularconeFit class constructor.
        %
        % Syntaxes:
        %
        %    (1) obj=rightcircularconeFit(xyz)
        %    (2) obj=rightcircularconeFit(xyz,v0)
        %                            
        %IN:                         
        %                            
        %      xyz: A 3xN  matrix whose columns are sample coordinates of a cylinder to be fitted. 
        %       v0: Initial guess of vertex location [vx0,vy0,vz0].
        %
        %NOTE: If syntax (1) is used, the routine will auto-generate an
        %      initial guess of the cylinder axis


            if ~nargin, XYZ=[]; end
            
            if isempty(XYZ), return; end
            
            
            XYZ0=XYZ;
            
            [xyzCtr,T]=quadricFit.homogNorm(XYZ);
            %T=eye(4); xyzCtr=XYZ; %DEBUGGUNG
            
            switch nargin
               
                case 1 %complete auto initialization
                    
   
                    warnStruct=warning('off');
                      obje=ellipsoidalFit(XYZ0);
                    warning(warnStruct);

                   v0=obje.center;
                    
                case 2 %inital vertex given     
                    
                    v0=varargin{1};
                
            end
            
            v0=v0(:);
            v0=rightcircularconeFit.tform(T,v0);
            
            vOpt=parameterEstimation(xyzCtr,v0);
            
            [~,dOpt,caOpt]=Fcost(vOpt,xyzCtr);
            
            
            R=quadricFit.vecrot([1;0;0],dOpt);
            vertex=quadricFit.invtform(T,vOpt);
            
            t=dOpt(:).'*(XYZ0-vertex);
            
            obj.XYZ=XYZ0;
            obj.vertex=vertex(:).';
            obj.cone_angle=caOpt;
            obj.ypr=quadricFit.rot2ypr(R);
            obj.height=max(t);
            
        end
        
        function xyz = sample(obj,varargin)
        %Generate discrete samples on fitted cylinder.
        %
        % xyz = obj.sample(azimSamps,axialSamps)
        %
        %IN:
        %
        %    azimSamps:  Azimuthal sample coordinates in degrees.
        %    axialSamps: Axial coordinates of samples (along along axis of
        %                cone with vertex at zero and positive elsewhere).
        %
        %OUT:
        %
        %    xyz: A 3x1 cell array of x,y,z sample locations
        %         on the fitted cylinder.
            
          assert(numel(varargin)<=2,'Too many input arguments');
          xyz = obj.xyzsim(obj.vertex,obj.cone_angle,obj.yp,varargin{:});
          xyz=num2cell(xyz,2);  
          
        end
        
        function varargout=fimplicit(obj,varargin)
            
            varargin=quadricFit.surfDefaults(varargin);
            
            x0=obj.vertex(1);
            y0=obj.vertex(2);
            z0=obj.vertex(3);

            ub=obj.height;
            lb=0;
            
            Rt=obj.R.';
            
            c=tand(obj.cone_angle);
            Q=diag([-c.^2,1,1]);
            
            hFit=fimplicit3(@quadform,varargin{:});
            
            if nargout, varargout={hFit}; end
            
            function d=quadform(x,y,z)
                
                xyz=Rt*[x(:).'-x0; y(:).'-y0; z(:).'-z0];
                
                inrange =  lb<=xyz(1,:) & xyz(1,:)<=ub;

                d=reshape( (sum((Q*xyz).*xyz)-1), size(x));
                
                d(~inrange)=nan;
            end
            
            
        end
        
        
        function varargout=surf(obj,varargin)
            
            varargin=quadricFit.surfDefaults(varargin);
            
            m=tand(obj.cone_angle);
            ub=obj.height;
            lb=0;
            
            
            delta=(ub-lb)/100;
            r=m*[ flip(delta:delta:abs(lb)), 0 ,delta:delta:abs(ub) ];
            r0=find(r==0,1);
            
            [X,Y,Z]=cylinder(r,1e3); 
            Z=Z*(ub-lb);
            Z=Z-Z(r0);
            
            T0=  [ 0     0     1
                  0     1     0
                 -1     0     0];
            
            T=[obj.R*T0,obj.vertex(:); 0 0 0 1];
            
            
            hFit=surf(X,Y,Z, 'Parent',hgtransform('Matrix',T),varargin{:});
            grid on
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
                propList = struct('vertex',obj.vertex,...
                                  'cone_angle', obj.cone_angle,...
                                  'height',obj.height,...
                                  'yaw',obj.yaw,...
                                  'pitch',obj.pitch);
                propgrp = matlab.mixin.util.PropertyGroup(propList);
            end
        end

    end
    
    methods (Static)
        
        function obj = groundtruth(xyz,vertex,cone_angle,height,yp)
        %Creata an object by specifying the geometry instead of fitting it.
        %Useful for plotting a ground truth surface.
        %
        % obj = rightcircularconeFit.groundtruth(xyz,vertex,cone_angle,height,yp)
        %
        %IN:
        %
        %       xyz:      a 3xN matrix of data points (can also be []) 
        %    vertex:      ground truth cone vertex [x0,y0,z0]
        %    cone_angle:  ground truth angle from cone axis to surface
        %    height:      ground truth cylinder height
        %        yp:  ground truth [yaw,pitch] of the cylinder in degrees
        
          
          obj=rightcircularconeFit([]);
          obj.XYZ=xyz;
          obj.vertex=vertex;
          obj.cone_angle=cone_angle;
          obj.height=height;
          obj.yp=yp;
                  
        end        
        
        
        
        function xyz = xyzsim(vertex,cone_angle,yp,theta,h,sig)
        % xyz = circularconeFit.xysim(vertex,cone_angle,yp,azimSamps,axSamps,sig)
        %
        %IN:
        %
        %    vertex:  vertex coordinates [vx,vy,vz]
        %    cone_angle:  cone angle in degrees (angle from cone axis
        %                 to surface).
        %        yp:  yaw, pitch angle vector [y,p] in degrees.
        %
        %    azimSamps:  Azimuthal sample coordinates in degrees.
        %    axialSamps: Axial coordinates of samples (along cone axis
        %                with zero at the vertex and positive elsewhere).
        %    sig: additive Gaussian noise standard deviation. 
        %
        %OUT:
        %
        %    xyz: A 3xN matrix whose columns are simulated sample locations
        %         on cylinder.
            
            
            if ~exist('sig','var')||isempty(sig), sig=0; end
            if ~exist('theta','var')||isempty(theta), theta=0:10:350; end
            if ~exist('h','var')||isempty(h), h=linspace(0,1,20); end
            
            assert(all(h>=0),'axialSamps must be non-negative')
            
            [theta,h]=ndgrid(theta,h);
            
            [X,Y,Z]=pol2cart(theta*pi/180,1,h);

            D=Z.*tand(cone_angle);
            
             xyz=[D(:).*X(:),D(:).*Y(:),Z(:)].';
            
            if sig
                
                N=size(xyz,2);
                %dr=normalize(xyz.*[1;1;0],1,'norm').*randn(1,N)*sig;
                dr=normalize(randn(3,N),1,'norm')*sig;
                
                xyz=xyz+dr; %add radial noise
                
            end
            
            xyz=xyz([3,1,2],:);
            
            
            R=quadricFit.ypr2rot([yp,0]);
            
            xyz=R*xyz+vertex(:);
            
        end
        
    end
end

%%%CLASS-RELATED FCNs

function [vOpt,fval]=parameterEstimation(xyzCtr,v0)



        fun=@(v) Fcost(v,xyzCtr);
        %opts=optimset('TolFun',1e-10,'MaxFunEvals',1e4, 'MaxIter',1e4);
        opts=optimset('TolFun',1e-10);
        

        [vOpt,fval]=fminsearch(fun, v0,opts);


end

function [fval,d,cone_angle]=Fcost(v,xyz)

      v=v(:);

      Zt=(xyz-v).';

      U=vecnorm(Zt,2,2);
      Unorm=norm(U);
      U=U/Unorm;

      A=Zt-U*(U.'*Zt);

      [~,S,V]=svd(A,0);

      fval=S(end);
      d=V(:,end);

      Ztd=Zt*d;
      if all(Ztd<=0)
          d=-d;
          Ztd=-Ztd;
      elseif ~all(Ztd>=0)
         fval=fval+expm1(-min(Ztd)); 
         %fval=inf;
      end

      if nargout<3, return; end

      c=(U.'*Ztd)/Unorm;

      cone_angle=acosd(c);


end
