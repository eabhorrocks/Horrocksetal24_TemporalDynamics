classdef linear2dFit<conicFit
%A class for executing and managing 2D line fits
    

   properties (SetAccess=protected)
       
       %angle; %inherited
       x1y1(2,1) double  %line segment end points
       x2y2(2,1) double
       
   end
    
    
    properties (Dependent,Hidden) 
        abc       (3,1)
        direction (1,2)
        normal    (1,2)
        p1,p2     (1,2)
        p1p2      (1,4)
    end
 
    properties (Dependent) %Not settable
        a  
        b  
        c
    end
    
    methods %Accessors
        
        function obj=set.p1p2(obj,p1p2)
           
            x1y1=p1p2(1:2).';
            x2y2=p1p2(3:4).';
            direction=x2y2-x1y1;
            direction=direction/norm(direction);
            obj.angle=atan2d(direction(2),direction(1));
            obj.x1y1=x1y1;
            obj.x2y2=x2y2;
            
        end
        
        
        function p1=get.p1(obj)
            p1=obj.x1y1(:).';
        end
        
        function p2=get.p2(obj)
            p2=obj.x2y2(:).';
        end
        
        function direction=get.direction(obj)
            
            direction=obj.R(:,1).';
            
        end
        
        
        function normal=get.normal(obj)
            
            normal=obj.R(:,2).';
            
        end
        
        function abc=get.abc(obj)
            
            el =cross([obj.x1y1;1],[obj.x2y2;1]);
            abc = el.'/norm(el(1:2));
             
        end
        
      
        
        function val=get.a(obj)
               val=obj.abc(1);
        end
        
        function val=get.b(obj)
               val=obj.abc(2);
        end
        
        function val=get.c(obj)
               val=obj.abc(3);
        end       
    end
    
    methods
        
        function obj=linear2dFit(XY)
        %linear2dFit constructor
        %
        %    obj=circularFit(xy)                    
        %                            
        %IN:                         
        %                            
        %    xy: a 2xN  matrix whose columns are sample coordinates on a line segment to be fitted.     
           
            if ~nargin, XY=[]; end

            if isempty(XY), return; end
            
            
            XY0=XY;
            
            [XY,T]=conicFit.homogNorm(XY);
            
            
            direction=conicFit.leastnull(XY.');
            
            if direction(1)<0 || (direction(1)==0 && direction(2)<0)
                direction = -direction;
            end
            

            t=direction.'*XY;

            x1y1=T\[min(t)*direction;1];
            x2y2=T\[max(t)*direction;1];
            
            obj.x1y1=x1y1(1:2)/x1y1(3);
            obj.x2y2=x2y2(1:2)/x2y2(3);
            obj.angle=atan2d(direction(2),direction(1));
            obj.XY=XY0;
 
        end
        

        
        
        function xy = sample(obj,t)
        %Generate discrete samples on the fitted 2D line.
        %
        %   xy = sample(obj,t)
        %
        %IN:      
        %
        %   t:  affine combination coefficient vector
        %
        %OUT:                        
        %                            
        %    xy: A 2x1 cell array of x,y sample locations
        %        along the line, given by(1-t(i))*obj.p1.'+t(i)*obj.p2.'
            
            
            xy=obj.xysim(obj.x1y1,obj.x2y2,t,0);
            xy=num2cell(xy,2);
          
        end
        

        function varargout=fimplicit(obj,varargin)
            
            varargin=conicFit.lineDefaults(varargin);
            
            abc=num2cell(obj.abc);
            [a,b,c]=deal(abc{:});
            pm=obj.x1y1/2+obj.x2y2/2;
            L=norm(pm-obj.x1y1);
            d=obj.direction;
            hFit=fimplicit( @fun , varargin{:});
            
            if nargout, varargout={hFit}; end

            function val=fun(x,y)
                
                p=abs(d*[x(:)'-pm(1);y(:)'-pm(2)])>L;
                val=(a*x+b*y+c);
                val(p)=nan;
                
            end
            
            
        end

        function varargout=line(obj,varargin)
        
             Args=num2cell( [obj.x1y1,obj.x2y2],2);
             
             hFit=line(Args{:},varargin{:});
             
             if nargout, varargout={hFit}; end
 
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
        
        function obj = groundtruth(xy,p1,p2)
        %Creata an object by specifying the geometry instead of fitting it.
        %Useful for plotting a ground truth locus.
        %
        % obj = linear2dFit.groundtruth(xy,p1,p2)
        %
        %IN:
        %
        %  xy: a 2xN matrix of data points (can also be []).
        %  p1:  one point on line [x1 y1]
        %  p2:  second point on line [x2 y2]   
     
            obj=linear2dFit([]);
            obj.XY=xy;
            obj.p1p2=[p1(:).',p2(:).'];
            
        end        
        
        function XY = xysim(p1,p2,t,sig)
        % xy = linear2dFit.xysim(p1,p2,t,sig)
        %
        %IN:                         
        %                            
        %  p1:  one point on line [x1 y1]
        %  p2:  second point on line [x2 y2]                 
        %
        %   t:  affine combination coefficeint vector
        %  sig: Gaussian noise sigma
        %
        %OUT:                        
        %                            
        %    xy: A 2xN matrix whose columns are simulated sample locations.
        %        With sig=0, they are given by(1-t(i))*p1.'+t(i)*p2.'
        
            
            
            if ~exist('sig','var')||isempty(sig), sig=0; end
            if ~exist('t','var')||isempty(t), t=linspace(0,1,20); end
            
            t=t(:).';
            
            XY=p1(:)*(1-t)+p2(:)*t;
            
            
            
            
            if sig
                
                N=numel(t);
                %B=null(p1-p2); dr=B*sig*randn(1,N);
                dr=normalize(randn(2,N),1,'norm')*sig;
                
                XY=XY+dr; %add radial noise
                
            end
            
            
        end
        
    end

    
end

