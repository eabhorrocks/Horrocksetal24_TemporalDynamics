classdef planarFit<quadricFit
%A class for executing and managing 3D plane fits
 
    properties (SetAccess=protected) 
        distance (1,1)
    end
 
    properties (Dependent) %Not settable
        normal   (1,3)
        a  
        b  
        c
        d %settable alias for distance
        closestPoint
    end
    
    methods %Accessors

        
        
        function normal=get.normal(obj)
            
            normal=obj.R(:,1).';
            
        end
        
        function obj=set.normal(obj,normal)
           
            normal=normal(:);
            s=sign(normal(find(normal,1,'last')));
            normal = s*normal/norm(normal);

            R=[normal,null(normal.')];

            R(:,3)=cross(R(:,1),R(:,2)); 
            obj.R=R;
            
        end

        function closestPoint = get.closestPoint(obj)
            
            closestPoint=obj.normal*obj.d;
            
        end
      
        
        function val=get.a(obj)
               val=obj.normal(1);
        end
        
        function val=get.b(obj)
               val=obj.normal(2);
        end
        
        function val=get.c(obj)
               val=obj.normal(3);
        end     
               
        function val=get.d(obj)
               val=obj.distance;
        end  
        
        function obj=set.d(obj,d)
               obj.distance=d;
        end       
    end
    
    methods
        
        function obj=planarFit(XYZ)
        %planarFit constructor
        %
        %    obj=planarFit(xyz)                    
        %                            
        %IN:                         
        %                            
        %    xyz: a 3xN  matrix whose columns are sample coordinates on a plane to be fitted.      
           
            if ~nargin, XYZ=[]; end

            if isempty(XYZ), return; end
            
            
            XYZ0=XYZ;
            
            [XYZ,T]=conicFit.homogNorm(XYZ);
            
            
            normal=conicFit.mostnull(XYZ.');
            normal=T(1:3,1:3)\normal;
             s=sign(normal(find(normal,1,'last')));
            normal = s*normal/norm(normal);
            
            d=dot([normal;0], T\[0;0;0;1]);
            
            R=[normal,null(normal.')];

            R(:,3)=cross(R(:,1),R(:,2));

            obj.XYZ=XYZ0;
            obj.R=R;
            obj.d=d;
 
        end
        

        
        
        function xyz = sample(obj,varargin)
        %Generate samples on the fitted plane.
        %
        %   xyz = obj.sample(t1,t2)
        %   xyz = obj.sample(b1,b2,t1,t2)
        %   xyz = obj.sample(b0,b1,b2,t1,t2)
        %
        %in:
        %
        %     b0: 3D position vector of a point on the plane.
        %  b1,b2: plane basis vectors. If either is set to [], it will be
        %         replaced by a vector orthogonal to the other.
        %  t1,t2: basis vector coefficients
        %
        %out:
        %
        %   xyz: A 3x1 cell array of x,y,z sample locations on fitted plane over 
        %        the parallelogram-shaped lattice b0+t1(i)*b1+t2(j)*b2 where b0 
        %        is the closest point on the plane to the 3D origin and b1 and b2
        %        are orthogonal plane basis vectors.
        %
        %Note: If b0, b1, and b2 do not span the fitted plane, they will be scaled
        %      and rotated so that they do so.
        

           N=numel(varargin);
           assert(N==2||N==4||N==5,'obj.sample(arg1,...,argN) must have N=2, N=4, or N=5 arguments');
           
           if N<5
            b0=obj.closestPoint;   
           else
            b0=varargin{1};  varargin(1)=[];
            b0=(  obj.d/dot( obj.normal(:),b0(:) ) )*b0 ; %scale so that b0 is in the plane
           end
           
           if N<3
              b1=[1,0,0];  b2=[];
           else
             [b1,b2]=deal(varargin{1:2});  varargin(1:2)=[]; 
           end           
        
           [t1,t2]=deal(varargin{:});
           
           if isempty(b1) && isempty(b2), 
               error 'At least one of b1 or b2 must be a non-empty 3-vector'; 
           elseif isempty(b2)
               b2=cross(obj.normal(:),b1(:));
           elseif isempty(b1)
               b1=cross(b2(:),obj.normal(:));
           end
 
           b1=b1(:);b2=b2(:);
           B=quadricFit.vecrot(cross(b1, b2) ,obj.normal(:))*[b1,b2] ;
           b1=B(:,1); 
           b2=B(:,2);

           xyz=obj.xyzsim(b0,b1,b2,t1,t2);
           xyz=num2cell(xyz,2);
            
        end
        

        function varargout=fimplicit(obj,varargin)
            
            varargin=quadricFit.surfDefaults(varargin);
            
            abc=num2cell(obj.normal);
            [a,b,c]=deal(abc{:});
            d=obj.d;
            
            hFit=fimplicit3( @(x,y,z) a*x+b*y+c*z - d , varargin{:});
            
            if nargout, varargout={hFit}; end
        end
        
        function varargout= showfit(obj,varargin)
            
            [varargout{1:nargout}]=obj.fimplicit(varargin{:});
            showfit@quadricFit(obj);
        end
        
    
    end
    
    
    methods (Access = protected, Hidden)
        function propgrp = getPropertyGroups(obj)
            if ~isscalar(obj)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
            else

                propList = struct(...
                    'normal',obj.normal,...
                    'distance',obj.d);
                propgrp = matlab.mixin.util.PropertyGroup(propList);
            end
        end
        

        
    end

    
    methods (Static)
        
        function obj = groundtruth(xyz,varargin)
        %Creata an object by specifying the geometry instead of fitting it.
        %Useful for plotting a ground truth surface.
        %
        %  (1) obj = planarFit.groundtruth(xyz,normal,d)
        %  (2) obj = planarFit.groundtruth(xyz,b0,b1,b2)
        %
        %IN:
        %
        % (1)
        %
        %         b0: 3D position vector of a point on the plane
        %      b1,b2: plane basis vectors
        %
        % (2)
        %
        %     normal: ground truth plane normal
        %          d: signed distance to the plane along the direction of
        %             the normal.
        
        
          
          obj=planarFit([]);
          obj.XYZ=xyz;

          
          switch nargin
              
              case 3 
                [normal,d]=deal(varargin{:});
                
              case 4
                
                  [b0,b1,b2]=deal(varargin{:});
                  normal=cross(b1(:),b2(:)); normal=normal/norm(normal);
                  d=dot(b0,normal);
                  
              otherwise
                  error 'Only 3 or 4 input arguments permitted'
                  
          end

          obj.normal=normal;
          obj.distance=d;
                  
        end        
        
        
        function XYZ = xyzsim(origin,basis1,basis2,t1,t2,sig)
        % xyz = planarFit.xyzsim(b0,b1,b2,t1,t2,sig)
        %
        %in:
        %
        %     b0: 3D position vector of a point on the plane
        %  b1,b2: plane basis vectors
        %  t1,t2: basis vector coefficients
        %  sig: Gaussian noise sigma
        %
        %out:
        %
        %   xyz: 3xN matrix of noisy sample locations. Before adding noise,
        %        the samples lie at locations b0+t1(i)*b1+t2(j)*b2 and
        %        therefore cover a paralleogram in the desired plane.
            
            if ~exist('sig','var')||isempty(sig), sig=0; end
            if ~exist('t1','var')||isempty(t1), t1=linspace(0,1,10); end
            if ~exist('t2','var')||isempty(t2), t2=linspace(0,1,10); end
            if ~exist('basis1','var')||isempty(basis1), basis1=[1,0,0]; end
            if ~exist('basis2','var')||isempty(basis2), basis2=[1,0,0]; end   
            if ~exist('origin','var')||isempty(origin), origin=[0,0,0]; end
            
            basis1=basis1(:); basis2=basis2(:); origin=origin(:);

            [T1,T2]=ndgrid(t1,t2);

            XYZ=origin+basis1*T1(:).'+ basis2*T2(:).';
            
            if sig
                
                normal=cross(basis1,basis2);
                normal=normal/norm(normal);
                
                N=size(XYZ,2);
                %dr=normal*(sig*randn(1,N));
                dr=normalize(randn(3,N),1,'norm')*sig;
                
                XYZ=XYZ+dr; %add radial noise
                
            end
            
            
        end
        
    end

    
end

