classdef Parameter
    
    properties
        % Problem parameters -------------------------------------------------
        pb_type
        % 101: Poission problem
        % 102: Lapalce Eigenvalue problem
        
        % 201: functional (u,g) with Poission problem
        % 202: Lapalce eigenvalue as a non-linear functional
        % 203?functional <q*n, phi>_T, where T is part of the boundary
        
        pb_parameters % all the necessary parameters to define the PDE,
                      % including the source_f, the boundary data 
        
        % domain and mesh parameters---------------------------------------
        
        structure_flag % 0: unstructured, 1: structured
        dom_type % 'Rec', 'L', 'Cir', 'CirHole'
        h0 
        dirichlet_flag % see different boundary names in Build2DMesh
        neuman_flag
        geo_parameters
        
        
        
        % numerical method parameters -------------------------------------
        
        order  % HDG method polynomial degree
        tau % stabiliaztion parameter
        post_process_flag
        
        % experiment parameters--------------------------------------------
        GQ_deg % Gauss Quadrature 
        
        Niter % how many mesh refinements
        
        refine_flag % 0: uniform refine, 1: adaptive refine
        
        report_flag % 1: report errory analysis if possibal
        
        visualize_flag % 1: visualize numerical solutions
        
        extra_parameters
        
        %------------------------------------------------------------------
    end
    
    methods 
        function obj = Parameter()
            cprintf('blue','Initiate empty Parameter object.\n')
        end
        
        function obj = SetPb(obj,pb_type,varargin)
            obj.pb_type = pb_type;
            obj.pb_parameters = varargin;  
            cprintf('blue','Set problem parameters... done\n')
        end
        
        function obj = SetMesh(obj,structure_flag,dom_type,h0,dirichlet_flag,neuman_flag,varargin)
            obj.structure_flag = structure_flag;
            obj.dom_type = dom_type;
            obj.h0 = h0;
            obj.dirichlet_flag = dirichlet_flag;
            obj.neuman_flag = neuman_flag;
            obj.geo_parameters = varargin;
            cprintf('blue','Set domain and mesh parameters... done\n')
        end
        
        function obj = SetNM(obj,order,tau, pp)
            obj.order = order;
            obj.tau = tau;
            obj.post_process_flag = pp;
            cprintf('blue','Set numerical methods parameters... done\n')
        end
        
        function obj = SetExp(obj, GQ_deg, Niter, refine_flag, report_flag,visualize_flag, varargin)
            obj.GQ_deg =GQ_deg;
            obj.Niter = Niter;
            obj.refine_flag = refine_flag;
            obj.report_flag = report_flag;
            obj.visualize_flag = visualize_flag;
            obj.extra_parameters = varargin;
            cprintf('blue','Set numericla experment parameters... done\n')
        end
        
        
        
    end
end
