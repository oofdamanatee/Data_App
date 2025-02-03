classdef movingBeamMdl < diffusion_function
    
    properties
        t
        gel_radius
        beam_size
        nmax
        rlim
        final_conc
        dx
        D
        t0
    end
    
    methods
        function obj = movingBeamMdl(options)
            obj@diffusion_function(options);
        end
        function obj = calcCurve(obj,p)
            obj = obj.updateFreeFitParams(p);
            uh = obj.paramStruct;
            y = diffusion_moving_beam(obj.t,uh.D,uh.gel_radius,...
                uh.final_conc,uh.nmax,uh.beam_size,uh.dx,0,'rlim',uh.rlim,...
                't0',uh.t0);
            obj.simMatrix = y';
        end
        function fh = getFunctionHandle(obj)
            str = "@(";
            for ii = 1:numel(obj.freeFitParamNames)
                str = str + obj.freeFitParamNames{ii} + ",";
            end
            str = str + "t)diffusion_moving_beam(t,obj.D.value,obj.gel_radius.value,obj.final_conc.value,obj.nmax.value,obj.sigma.value,obj.dx.value,0,""rlim"",obj.rlim.value)";
            fh = str2func(str);
        end
    end
    
end