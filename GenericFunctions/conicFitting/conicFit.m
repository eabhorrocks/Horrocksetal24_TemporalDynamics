classdef (Abstract) conicFit < matlab.mixin.CustomDisplay
%A base class for executing and managing 2D conic fits
    
    properties (SetAccess=protected)
        XY
        angle=0; 
    end
    properties (Dependent,Hidden)
        R
        theta;
    end
    
    
    methods

        function R=get.R(obj)
            x=obj.angle;
            R=[cosd(x),-sind(x); sind(x),cosd(x)];
        end
        
        function obj=set.theta(obj,val)
            obj.angle=val;
        end
        function val=get.theta(obj)
            val=obj.angle;
        end
        
        function varargout=scatter(obj,varargin)
        %Scatter plot of the sample coordinate data.
        %
        %  hData=scatter(obj,Name,Value)
        %
        %The Name/Value pairs are passed as input to scatter()
        
            if isempty(obj.XY), [varargout{1:nargout}]=deal([]);return; end   
        
            varargin= [{'filled','SizeData',50},varargin]; 
            
            %varargout{1}=scatter(obj.XY(1,:),obj.XY(2,:),varargin{:});
            h=scatter(obj.XY(1,:),obj.XY(2,:),varargin{:});
            
             xlabel 'X'; ylabel 'Y';
            
            if nargout, varargout={h};end
            
        end
        
        function varargout=plot(obj,varargin)
        %Plot the resulting fit (fimplicit) overlaid with scatter plot of 
        %sample points.
        %
        %   [hFit,hData]=plot(obj,fitArgs={},scatArgs={})
        %                           
        %OUT:                        
        %                            
        %     hFit: handle to the fitted curve fimplicit plot                 
        %    hData: handle to scatter plot of samples                  
        %                            
        %IN:                         
        %
        %   fitArgs: cell array of Name/Value options for fimplicit plot
        %   scatArgs: cell array of Name/Value options for scatter plot.              
   


            
            funfit=@obj.showfit;
            [varargout{1:nargout}]=obj.verify(funfit,varargin{:});

        end 
        
        function varargout=verify(obj,funfit, fitArgs, scatArgs)
        %Plot the resulting fit overlaid with scatter plot of 
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
        %     funfit: function handle to a command for plotting the fit.              
        %    fitArgs: cell array of Name/Value options for funfit().
        %   scatArgs: cell array of Name/Value options for scatter plot.
        
            if nargin<4,scatArgs={}; end
            if nargin<3,fitArgs={}; end
            
            hData=obj.scatter(scatArgs{:});
            axis equal
            
            hold on;
            hFit=funfit(fitArgs{:});
            hold off 
            
            if nargout, varargout={hFit,hData}; end
            
        end
        
        function varargout=showfit(obj,varargin)
              [varargout{1:nargout}]=obj.fimplicit(varargin{:});
              xlabel X; ylabel Y;
        end
        
    end
    
    methods (Abstract)
        
       varargout = fimplicit(obj,varargin);
        
    end
    methods  (Static, Abstract)
        
       obj  = groundtruth(varargin); 
       hFit = xysim(obj,varargin);
        
    end    
    
    methods (Static,Hidden)
        function [zn,T]=homogNorm(z)
            %Utility function for normalizing/centering homogeneous coordinates
            
            zc=mean(z,2);
            
            d=length(zc); % dimension of z(:,i)
            
            
            zn=bsxfun(@minus,z,zc);
            
            s=sqrt(d)./mean(vecnorm(zn,2,1));
            
            zn=zn*s;
            
            A=diag([s*ones(d,1);1]);
            B=eye(size(A)); B(1:end-1,end)=-zc(:);
            T=A*B;
        end
        
         function x=mostnull(A)
            
            [~,~,V]=svd(A,0);
            x=V(:,end);
            
        end
        
        function x=leastnull(A)
            [~,~,V]=svd(A,0);
            x=V(:,1);
        end 
        
        function dcell=lineDefaults(dcell)
            
             dcell=[{'LineWidth',1.5,'Color','r'},dcell];
            
        end
        
       function zn=tform(T,z)
            
            zn=T(:,1:end-1)*z+T(:,end);
             zn(zn(end,:)==0)=0;
            zn=zn(1:end-1,:)./zn(end,:); 
            
        end
         function zn=invtform(T,z)

            zn=conicFit.tform(inv(T),z);
            
        end      

    end
end

