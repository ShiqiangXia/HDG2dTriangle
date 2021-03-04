classdef Mesh
    properties
        dom_type; % domain type 'square', 'Lshape'
        mesh_type; % triangular or quad
        
        num_elements;
        num_faces;
        num_verices;
        
        element_list;  %Num_elements x 3(tri)/4(quad), each row is the index of the vertices of this element
        face_list; %Num_faces x 2, each row is the index of the vertices of the face(edge)
        vertices_list; %Num_verices x 2, coordinates of vertices
        
        element_faces_list; % the index of the faces of the element
        
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
            
            figure;
            
            p = obj.vertices_list;
            e = obj.element_list;
            faces_dirichlet = find(obj.f_type==1);
            faces_neuman = find(obj.f_type==2);
            
            nodes_dirichlet = [obj.face_list(faces_dirichlet,1);obj.face_list(faces_dirichlet,2)];
            nodes_dirichlet = unique(nodes_dirichlet);
            
            nodes_neuman = [obj.face_list(faces_neuman,1);obj.face_list(faces_neuman,2)];
            nodes_neuman = unique(nodes_neuman);
            
            nvertice = length(p);
            labs = 1:nvertice;
            triplot(e,p(:,1),p(:,2));
            hold on;
            plot(p(:,1),p(:,2),'k.',...
                p(nodes_dirichlet,1),p(nodes_dirichlet,2),'rx',...
                p(nodes_neuman,1),p(nodes_neuman,2),'bd');
    
            legend('faces','vertices','dirichlet', 'neuman');
            
            if nvertice < 150 && vertex_flag
                labelpoints(p(:,1),p(:,2),labs,'NW',0.5,0,'FontSize', 14);
            end
            %simpplot(obj.vertices_list,obj.element_list);
        end
        
        function mymesh = Refine(obj, marked, flag)
            % step 1: Refine the marked elements
            p = obj.vertices_list;
            e = obj.element_list;
            if strcmp(flag, 'RGB')
                [new_p,new_e] =TrefineRGB(p,e,marked); 
            elseif strcmp(flag, 'RG')
                [new_p,new_e] =TrefineRG(p,e,marked);
                
            elseif strcmp(flag, 'R')
                fprintf("Warning: R refinement makes hanging nodes!\n")
                [new_p,new_e] =TrefineR(p,e,zeros(0,3),marked);
            elseif strcmp(flag, 'NVB')
                [new_p,new_e] =TrefineNVB(p,e,marked);
            else
                error('Refine flag is not correct!');
            end
            
            new_e = Countclockwise(new_p,new_e); % make sure a counterclockwise ordering of the vertices

            [f,ef,face_type] = LabelFaces(obj.dom_type,new_p,new_e,obj.dirichlet_flag,obj.neuman_flag,obj.dom_parameters{:});
           
            mymesh = Mesh(obj.dom_type,new_p,new_e,f,ef,face_type,obj.dirichlet_flag,obj.neuman_flag,obj.dom_parameters{:});
                
        end
        
        
    end
    
end