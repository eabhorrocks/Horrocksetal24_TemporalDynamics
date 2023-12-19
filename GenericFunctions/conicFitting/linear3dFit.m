classdef linear3dFit < matlab.mixin.CustomDisplay
%A class for executing and managing 3D line fits
    

   properties (SetAccess=protected)
       XYZ
       x1y1z1(3,1) double  %line segment end points
       x2y2z2(3,1) double
       
   end
    
    
    properties (Dependent,Hidden)
        p1,p2
        p1p2
        direction
    end
 
    
    methods %Accessors
        
        
        function obj=set.p1p2(obj,p1p2)
           
            x1y1z1=p1p2(1:3).';
            x2y2z2=p1p2(4:6).';

            obj.x1y1z1=x1y1z1;
            obj.x2y2z2=x2y2z2;
            
        end
        
        
        function direction=get.direction(obj)
            
            direction=(obj.x2y2z2-obj.x1y1z1).';
            
        end
        
        function p1=get.p1(obj)
              p1=obj.x1y1z1(:).';
        end
        
        function p2=get.p2(obj)
              p2=obj.x2y2z2(:).';
        end
        
    end
    
    methods
        
        function obj=linear3dFit(XYZ)
        %linear3dFit constructor
        %
        %    obj=linear3dFit(xyz)                    
        %                            
        %IN:                         
        %                            
        %    xyz: a 3xN  matrix whose columns are sample coordinates on a line to be fitted.     
           
            if ~nargin, XYZ=[]; end

            if isempty(XYZ), return; end
            
            
            XYZ0=XYZ;
            
            [XYZ,T]=conicFit.homogNorm(XYZ);
            
            
            direction=quadricFit.leastnull(XYZ.');
            
            idx=find(direction,1);
            if direction(idx)<0
                direction = -direction;
            end
            

            t=direction.'*XYZ;

            x1y1z1=T\[min(t)*direction;1];
            x2y2z2=T\[max(t)*direction;1];
            
            obj.x1y1z1=x1y1z1(1:3)/x1y1z1(4);
            obj.x2y2z2=x2y2z2(1:3)/x2y2z2(4);
            obj.XYZ=XYZ0;
 
        end
        
        function varargout=scatter(obj,varargin)
        %Scatter plot of the sample coordinate data.
        %
        %  hData=scatter(obj,Name,Value)
        %
        %The Name/Value pairs are passed on as input to scatter3()
        
            if isempty(obj.XYZ), [varargout{1:nargout}]=deal([]);return; end   
        
            varargin= [{'filled','SizeData',50},varargin]; 
            
            h=scatter3(obj.XYZ(1,:),obj.XYZ(2,:),obj.XYZ(3,:),varargin{:});
            
            [varargout{1:nargout}]=deal(h);
            
            xlabel X; ylabel Y; zlabel Z;
        end
        
        function varargout=plot(obj,varargin)
        %Plot the resulting fit (fimplict) overlaid with scatter plot of 
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
        %   fitArgs: cell array of Name/Value options for fimplicit plot
        %   scatArgs: cell array of Name/Value options for scatter3 plot. 
            

            
            funfit=@obj.fimplicit;
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
            if nargout, varargout={hFit,hData}; end
            
        end
              
        
        function xyz = sample(obj,t)
        %Generate samples on the fitted 3D line.
        %
        %    xy = obj.sample(t)
        %
        %IN:
        %
        %   t: affine combination coefficient values
        %
        %OUT:
        %
        %    xyz: 3x1 cell array of x,y,z sample locations on the fitted 3D
        %         line, given by(1-t(i))*obj.p1(:) + t(i)*p2(:)
            
            
            xyz=obj.xysim(obj.p11,obj.p2,t);
            xyz=num2cell(xyz,2);
            
            
        end
        
        function varargout=fimplicit(obj,varargin)
        
            xmym=(obj.x2y2z2+obj.x1y1z1)/2;
            
            L=(obj.x2y2z2-obj.x1y1z1)/2;
            
            [~,theta,ax]=quadricFit.vecrot([norm(L),0,0],L);
            
            T=makehgtform('translate',xmym,'axisrotate',ax,theta*pi/180);
            
            hFit=fimplicit(@fun,'Parent',hgtransform('Matrix',T) , varargin{:});
            
            if nargout, varargout={hFit}; end
            xlabel X; ylabel Y; zlabel Z;
            function val=fun(x,y)
                
                p=abs(x)>norm(L);
                val=y;
                val(p)=nan;
                
            end
            
        end
        
         function varargout=line(obj,varargin)
        
             Args=num2cell( [obj.x1y1z1,obj.x2y2z2],2);
             
             hFit=line(Args{:},varargin{:});
             
             if nargout, varargout={hFit}; end
             xlabel X; ylabel Y; zlabel Z;
         end
         
        function varargout= showfit(obj,varargin)
            
            [varargout{1:nargout}]=obj.line(varargin{:});
            xlabel X; ylabel Y; zlabel Z;
        end   
         
    end
    
    
    methods (Access = protected, Hidden)
        function propgrp = getPropertyGroups(obj)
            if ~isscalar(obj)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
            else
                propList = struct(...
                    'p1',obj.p1,...
                    'p2',obj.p2);
                propgrp = matlab.mixin.util.PropertyGroup(propList);
            end
        end
        

        
        
        
    end

    
    methods (Static)
        
        function obj = groundtruth(xyz,p1,p2)
        %Creata an object by specifying the geometry instead of fitting it.
        %Useful for plotting a ground truth locus.
        %
        % obj = linear3dFit.groundtruth(xy,p1,p2)
        %
        %IN:
        %
        %  xyz: a 3xN matrix of data points (can also be []).
        %  p1:  one point on line [x1 y1 z1]
        %  p2:  second point on line [x2 y2 z2]   
     
            obj=linear3dFit([]);
            obj.XYZ=xyz;
            obj.p1p2=[p1(:).',p2(:).'];
            
        end              
        
        
        
        
        function XYZ = xyzsim(p1,p2,t,sig)
        % xy = linear3dFit.xysim(p1,p2,t,sig)
        %
        %in:
        %
        %  p1: one point on line [x1 y1 z1]
        %  p2: second point on line [x2 y2 z2]
        %   t: affine combination coefficient values
        %  sig: Gaussian noise sigma
        %
        %out:
        %
        %    xyz: A 3xN matrix whose columns are simulated sample locations.
        %         With sig=0, they are given by(1-t(i))*p1.'+t(i)*p2.'
            
            
            if ~exist('sig','var')||isempty(sig), sig=0; end
            if ~exist('t','var')||isempty(t), t=linspace(0,1,20); end
            
            t=t(:).';
            
            XYZ=p1(:)*(1-t)+p2(:)*t;
            
            
            
            
            if sig
                
                N=numel(t);
                %B=null(p1-p2); dr=B*sig*normalize( randn(2,N) ,1, 'norm');
                dr=normalize(randn(3,N),1,'norm')*sig;
                
                XYZ=XYZ+dr; %add radial noise
                
            end
            
            
        end
        
    end

    
end

