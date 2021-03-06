classdef Mesh
    properties
        dom_type; % domain type 'Rec', 'L', 'Cir', 'CirHole'
        mesh_type; % triangular or quad
        
        num_elements;
        num_faces;
        num_verices;
        
        element_list;  %Num_elements x 3(tri)/4(quad), each row is the index of the vertices of this element
        face_list; %Num_faces x 2, each row is the index of the vertices of the face(edge)
        vertices_list; %Num_verices x 2, coordinates of vertices
        
        element_faces_list; % the index of the faces of the element
        
        Jacobian_list % the jacobian for mapping to the ref element 
        
        uhat_dir_list
        
        f_type; % Num_facesx1, record the type of each face
        
        % 0: interior
        % 1: Dirichlet
        % 2: Neuman
        
       
        dirichlet_flag;
        neuman_flag;
        
        dom_parameters;
        
        

    end
    
    methods
        function obj= Mesh(type,p,e,f,ef,f_type,dirichlet_flag,neuman_flag,varargin)
            
                
            obj.dom_type = type;
            
            obj.num_verices = length(p);
            [ne,nf] = size(e);
            
            obj.num_elements = ne;
            
            obj.Jacobian_list = TriJacobian(p,e);
            
            obj.uhat_dir_list = GetUhatOritation(e);
            
            if nf == 3
                obj.mesh_type = 'Tri';
            elseif nf == 4
                obj.mesh_type = 'Quad';
            else
                error('wrong mesh type.')
            end
            
            obj.num_faces = length(f);
            
            obj.element_list = e;
            obj.face_list = f;
            obj.vertices_list = p;
            
            obj.element_faces_list = ef;
            
            obj.f_type = f_type;
            
            obj.dirichlet_flag = dirichlet_flag;
            obj.neuman_flag = neuman_flag;
            if nargin - 8 > 0
                
                obj.dom_parameters = varargin;
            else
                
                obj.dom_parameters = {};
            end
            
            if 0
                obj.Plot(1);
            end
             
            
            
        end
        
        function Plot(obj,vertex_flag)
            
            bcol=[.8,.9,1];
            
            figure;
            
            p = obj.vertices_list;
            e = obj.element_list;
            
            face_type = unique(obj.f_type);
            % face_type(1) is inner face
            face_type1 = find(obj.f_type==face_type(2));
            
            
            nodes_type1 = [obj.face_list(face_type1,1);obj.face_list(face_type1,2)];
            nodes_type1 = unique(nodes_type1);
            
            if length(face_type)>2
                face_type2 = find(obj.f_type==face_type(3));
            else
                face_type2 = [];
            end
            nodes_type2 = [obj.face_list(face_type2,1);obj.face_list(face_type2,2)];
            nodes_type2 = unique(nodes_type2);
            
            nvertice = length(p);
            labs = 1:nvertice;
            trimesh(e,p(:,1),p(:,2),0*p(:,1),'facecolor',bcol,'edgecolor','k');
            hold on;
            plot(p(:,1),p(:,2),'k.')
            legend('element','vertices');
            
            % mark the dirichlet and neuman boundary nodes
            plot(p(:,1),p(:,2),'k.',...
                p(nodes_type1,1),p(nodes_type1,2),'rx',...
                p(nodes_type2,1),p(nodes_type2,2),'bd')     
            if isempty(nodes_type1) && isempty(nodes_type2)
                legend('element','vertices');
            elseif isempty(nodes_type1)
                legend('element','vertices', 'neuman')
            elseif isempty(nodes_type2)
                legend('element','vertices', 'dirichlet')
            else
                legend('element','vertices','bdry 1', 'bdry 2');
            end
           
            
            if   vertex_flag && nvertice < 150
                labelpoints(p(:,1),p(:,2),labs,'NW',0.5,0,'FontSize', 14);
            end
            view(2)
            axis equal
            ax=axis;axis(ax*1.001);
            if vertex_flag && nvertice>=150
                cprintf('UnterminatedStrings', '%d vertices are too many to plot, so I ignored them.\n',nvertice);
            end
            %simpplot(obj.vertices_list,obj.element_list);
            
            
        end
        
        function mymesh = Refine(obj, marked, flag)
            % step 1: Refine the marked elements
            p = obj.vertices_list;
            e = obj.element_list;
            face_type = unique(obj.f_type);
            
            face_id1 = face_type(2);
            
            if length(face_type)>2
                face_id2 = face_type(3);
            else
                face_id2 = 2;
            end
            faces_dirichlet = find(obj.f_type==face_id1);
            faces_neuman = find(obj.f_type==face_id2);
            
            diri = [obj.face_list(faces_dirichlet,1), obj.face_list(faces_dirichlet,2)];
            neu = [obj.face_list(faces_neuman,1), obj.face_list(faces_neuman,2)];
            
            if strcmp(flag, 'RGB')
                [new_p,new_e,diri,neu] =TrefineRGB(p,e,diri,neu,marked); 
            elseif strcmp(flag, 'RG')
                [new_p,new_e,diri,neu] =TrefineRG(p,e,diri,neu,marked);
                
            elseif strcmp(flag, 'R')
                cprintf('*UnterminatedStrings','Warning: R refinement makes hanging nodes!\n')
                [new_p,new_e,diri,neu] =TrefineR(p,e,zeros(0,3),diri,neu,marked);
                
            elseif strcmp(flag, 'NVB')
                [new_p,new_e,diri,neu] =TrefineNVB(p,e,diri,neu,marked);
            else
                error('Refine flag is not correct!');
            end
            
            dirinodes = unique(reshape(diri,[],1));
            neunodes = unique(reshape(neu,[],1));
            
            new_e = Countclockwise(new_p,new_e); % make sure a counterclockwise ordering of the vertices
            new_e = LongestEdgeFirst(new_p,new_e);% make sure the longest edge is the first element 
            % (this is used for mesh refinement to keep the shape-regularity)
            
            bdry_edges = boundedges(new_p,new_e);
            [f,ef,face_type] = LabelFaces(new_e,bdry_edges, dirinodes,face_id1,neunodes,face_id2);
            if ~(strcmp(obj.dom_type,'Rec') || strcmp(obj.dom_type,'L') || strcmp(obj.dom_type,'Comp_Rec'))
                cprintf('*UnterminatedStrings', 'Warning:\n')
                cprintf('UnterminatedStrings','Mesh refinment on a curved domain (%s) will NOT guarantee a better resolution of the curve!\n',obj.dom_type);
                cprintf('UnterminatedStrings','You many consider using Build2DMesh to generate a finer mesh to resolve the curved bounaries.\n')
            end
            mymesh = Mesh(obj.dom_type,new_p,new_e,f,ef,face_type,obj.dirichlet_flag,obj.neuman_flag,obj.dom_parameters{:});
                
        end
        
        function mymesh = UniformRefine(obj)
            n = obj.num_elements;
            marked = 1:n;
            mymesh = obj.Refine(marked,'RGB');
 
        end
        
        function PlotElement(obj,label_flag)
            bcol=[.8,.9,1];
            
            figure;
            
            p = obj.vertices_list;
            e = obj.element_list;
          
            
            % plot the mesh
            trimesh(e,p(:,1),p(:,2),0*p(:,1),'facecolor',bcol,'edgecolor','k');
            hold on;
            plot(p(:,1),p(:,2),'k.')
            legend('element','vertices');
            
            num = obj.num_elements;
            labs = 1:num;
            e1 = e(:,1);
            e2 = e(:,2);
            e3 = e(:,3);
            
            ave_x = (p(e1,1) + p(e2,1) + p(e3,1))/3.0;
            ave_y = (p(e1,2) + p(e2,2) + p(e3,2))/3.0;
            
          
            
            if  num < 250 && label_flag == 1
                labelpoints(ave_x,ave_y,labs,'center','FontSize', 14);
            end
            
            view(2)
            axis equal
            ax=axis;axis(ax*1.001);
            if num >= 250 && label_flag == 1
                labelpoints(ave_x(1:250),ave_y(1:250),labs(1:250),'center','FontSize', 14);
                cprintf('UnterminatedStrings', '%d elements are too many to label, so I only plot the first 250.\n',num);
            end
            %simpplot(obj.vertices_list,obj.element_list);
            
        end
        
        function Plot2(obj,vertex_flag, title_text)
            bcol=[.8,.9,1];
            
            figure;
            
            p = obj.vertices_list;
            e = obj.element_list;
           
            nvertice = length(p);
            labs = 1:nvertice;
            trimesh(e,p(:,1),p(:,2),0*p(:,1),'facecolor',bcol,'edgecolor','k');
            hold on;
            plot(p(:,1),p(:,2),'k.')
            legend('element','vertices');
            title(title_text)
           
            
            if   vertex_flag && nvertice < 150
                labelpoints(p(:,1),p(:,2),labs,'NW',0.5,0,'FontSize', 14);
            end
            view(2)
            axis equal
            ax=axis;axis(ax*1.001);
            
            if vertex_flag && nvertice>=150
                cprintf('UnterminatedStrings', '%d vertices are too many to plot, so I ignored them.\n',nvertice);
            end
            %simpplot(obj.vertices_list,obj.element_list);
            
        end
    end
    
end