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
        
        BCType;
        dirichlet_bdry % face index of dirichlet boundary
        neuman_bdry % face index of neuman boundary
    end
    
    methods
        function obj=Mesh(type,p,e,f,ef,dirichlet,neuman,BCType)
            if nargin>0
                
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
            
            obj.BCType = BCType;
            obj.dirichlet_bdry = dirichlet;
            obj.neuman_bdry = neuman;
            
            end
            
        end
        
        function mesh_plot(obj)
            simpplot(obj.vertices_list,obj.element_list);
        end
        
        
    end
    
end