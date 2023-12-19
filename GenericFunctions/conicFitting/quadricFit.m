classdef (Abstract) quadricFit < matlab.mixin.CustomDisplay
%A base class for executing and managing 3D quadric fits
    
    properties (SetAccess=protected)
        XYZ
        yaw=0;
        pitch=0;
        roll=0;
    end

    properties (Dependent,Hidden)
        R
        y,p,r
        ypr
    end
    
 
    methods
        
        function R=get.R(obj)
            
            R=quadricFit.ypr2rot(obj.ypr);
            
        end
        
        function obj=set.R(obj,R)
            
            
            obj.ypr=quadricFit.rot2ypr(R);
            
        end    
        
        function obj=set.ypr(obj,ypr)
           
            obj.yaw=ypr(1); obj.pitch=ypr(2); obj.roll=ypr(3);
            
        end
        
        function y=get.y(obj)
           y=obj.yaw;
        end
        function p=get.p(obj)
           p=obj.pitch;
        end       
        function r=get.r(obj)
           r=obj.roll;
        end        
        function ypr=get.ypr(obj)
             ypr=[obj.yaw,obj.pitch,obj.roll]; 
        end
        
        function obj=set.y(obj,y)
           obj.yaw=y;
        end
        function obj=set.p(obj,p)
           obj.pitch=p;
        end       
        function obj=set.r(obj,r)
           obj.roll=r;
        end 


        function varargout=scatter(obj,varargin)
        %Scatter plot of the sample coordinate data.
        %
        %  hData=scatter(obj,Name,Value)
        %
        %The Name/Value pairs are passed on as input to scatter3()
        
            if isempty(obj.XYZ), [varargout{1:nargout}]=deal([]);return; end   
        
            varargin= [{'filled','SizeData',50},varargin]; 
            
            %varargout{1}=scatter3(obj.XYZ(1,:),obj.XYZ(2,:),obj.XYZ(3,:),varargin{:});
            h=scatter3(obj.XYZ(1,:),obj.XYZ(2,:),obj.XYZ(3,:),varargin{:});
            
            xlabel 'X'; ylabel 'Y', zlabel 'Z'
            
            if nargout, varargout={h};end
        end
        
        function varargout=plot(obj,varargin)
        %Plot the resulting fit (with fimplicit3) overlaid with scatter3 plot of 
        %sample points.
        %
        %   [hFit,hData]=plot(obj,scatArgs={},fitArgs={})
        %
        %                            
        %OUT:                        
        %                            
        %     hFit: handle to the fitted curve fimplicit plot                 
        %    hData: handle to scatter plot of samples                  
        %                            
        %IN:                         
        %
        %   fitArgs: cell array of Name/Value options for fimplicit3 plot
        %   scatArgs: cell array of Name/Value options for scatter3 plot. 
            

            
            funfit=@obj.showfit;
            [varargout{1:nargout}]=obj.verify(funfit,varargin{:});
             xlabel 'X'; ylabel 'Y', zlabel 'Z'
        end 
        
        function varargout=verify(obj, funfit, fitArgs, scatArgs)
        %Plot the resulting fit overlaid with scatter3 plot of 
        %sample points.
        %
        %   [hFit,hData]=verify(obj, funfit,scatArgs={},fitArgs={})
        %
        %                            
        %OUT:                        
        %           
        %     hFit:   handle to the fitted curve plot                 
        %    hData:   handle to scatter plot of samples                  
        %                            
        %IN:                         
        %
        %     funfit: function handle to a command for plotting the fit
        %    fitArgs: cell array of Name/Value options for funfit().
        %   scatArgs: cell array of Name/Value options for scatter3 plot.  
            
            
            if nargin<3,fitArgs={}; end
            if nargin<4,scatArgs={}; end
 
            
            hData=obj.scatter(scatArgs{:});
            axis equal
            
            hold on;
            hFit=funfit(fitArgs{:});
            hold off 
            
            view(3)
            grid on
            
            xlabel 'X'; ylabel 'Y', zlabel 'Z'
             
            if nargout, varargout={hFit,hData}; end
            
        end
        
        function showfit(~) 
            xlabel X; ylabel Y; zlabel Z;
        end
        
    end
    
    methods (Abstract)
        
       varargout = fimplicit(obj,varargin);


       
    end

    
    methods  (Static, Abstract)
        
       obj  = groundtruth(varargin); 
       hFit = xyzsim(obj,varargin);
        
    end   
    
    
    methods (Static, Hidden)
        function [zn,T]=homogNorm(z)
            %Utility function for normalizing/centering homogeneous coordinates
            %
            %  [zn,T]=homogNorm(z)
            %
            % zn - normalized columns of z
            %  T - forward normalization transform
            
            zc=mean(z,2);
            
            d=length(zc); % dimension of z(:,i)
            
            
            zn=bsxfun(@minus,z,zc);
            
            s=sqrt(d)./mean(vecnorm(zn,2,1));
            
            zn=zn*s;
            
            A=diag([s*ones(d,1);1]);
            B=eye(size(A)); B(1:end-1,end)=-zc(:);
            T=A*B;
            
        end
        
        function ypr=rot2ypr(R0)
        %Map a 3D rotation matrix to yaw-pitch-roll (intrinsic z-y'-x'' rotation)
        %
        % [ypr,quat]=rot2ypr(R)
        %
        %in:
        %
        % R: A 3D rotation matrix
        %
        %out:
        %
        % ypr: A vector of yaw-pitch-roll angles in degrees [yaw,pitch,roll]
            
            R=R0;
            
            u=R(:,1);
            u=u/norm(u);
            [az,el]=cart2sph(u(1),u(2),u(3));
            y=az*180/pi; p=-el*180/pi;
            
            Ryp=quadricFit.ypr2rot([y,p,0]);
            
            Rr=R*Ryp.';
            
            
            Rskew=(Rr-Rr.')/2;
            v=Rskew([6,3,2]); v(2)=-v(2);
            ctheta=trace(Rr)/2-0.5;
            stheta=v*u;
            
            ypr=[y,p,atan2d(stheta,ctheta)];
            
            
            
        end
        
        function R=ypr2rot(ypr)
            
            R=eye(3);
            
            for i=1:3
                
                R= quadricFit.R3d(ypr(i), R(:,4-i))*R;
                
            end
            
            
        end
        
        function [R,theta,Nrm]=vecrot(vstart,vend)
        %Find rotation carrying one vector toward another about their common perpendicular
        %axis.
        %
        %IN:
        %
        %    vstart: Initial vector
        %      vend: Final vector
        %
        %OUT:
        %
        %    R: rotation matrix carrying vstart to vend in a rotation about the axis Nrm=cross(vstart,vend)
        %    theta:  the rotation angle in degrees
        %    Nrm: the rotation axis
            
            
            
            vstart=vstart(:)/norm(vstart);
            vend=vend(:)/norm(vend);
            
            Nrm=cross(vstart,vend);
            
            b=vend.'*vstart;
            
            theta = atan2d(sqrt(max(1-b^2,0)),b);
            
            R=quadricFit.R3d(theta,Nrm);
        end
        
        function R=R3d(deg,u)
        %R3D - 3D Rotation matrix counter-clockwise about an axis.
        %
        %R=R3d(deg,axis)
        %
        %  deg:    The counter-clockwise rotation about the axis in degrees.

           if deg==0, R=eye(3); return; end
            
            R=eye(3);
            u=u(:)/norm(u);
            x=deg; %abbreviation
            
            for ii=1:3
                
                v=R(:,ii);
                
                R(:,ii)=v*cosd(x) + cross(u,v)*sind(x) + (u.'*v)*(1-cosd(x))*u;
                %Rodrigues' formula
                
            end
            
        end
        
        
        function x=mostnull(A)
            
            [~,~,V]=svd(A,0);
            x=V(:,end);
            
        end
        
        function x=leastnull(A)
            [~,~,V]=svd(A,0);
            x=V(:,1);
        end
        
        
        function dcell=surfDefaults(dcell)
            
             dcell=[{'EdgeColor','none','FaceColor','r','FaceAlpha',0.3},dcell];
            
        end
        
        function zn=tform(T,z)
            
            zn=T(:,1:end-1)*z+T(:,end);
             zn(zn(end,:)==0)=0;
            zn=zn(1:end-1,:)./zn(end,:); 
            
        end
         function zn=invtform(T,z)

            zn=quadricFit.tform(inv(T),z);
            
        end       
        
    end
end

%%%%CLASS-RELATED FUNCS


