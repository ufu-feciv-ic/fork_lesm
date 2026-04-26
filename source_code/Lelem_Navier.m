%% Lelem_Navier Class
%
% This is a sub-class, in the Object Oriented Programming (OOP) paradigm,
% of super-class <lelem.html *Lelem*> in the <main.html LESM> (Linear Elements
% Structure Model) program. This sub-class implements abstract methods,
% declared in super-class *Lelem*, that deal with Navier (Euler-Bernoulli)
% flexural behavior of linear elements subjected to loading effects.
%
%% Navier (Euler-Bernoulli) Beam Theory
%
% In Euler-Bernoulli flexural behavior, it is assumed that there is no
% shear deformation. As consequence, bending of a linear element is such that
% its cross-section remains plane and normal to the element longitudinal axis.
%
%% Class definition
% Definition of sub-class *Lelem_Navier* derived from super-class <lelem.html Lelem>.
classdef Lelem_Navier < Lelem
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function load = Lelem_Navier(elem,unifdir,unifgbl,uniflcl,lindir,lingbl,linlcl,tvx,tvy,tvz)
            load = load@Lelem(elem);
            
            if (nargin > 1)
                load.uniformDir = unifdir;
                load.uniformGbl = unifgbl;
                load.uniformLcl = uniflcl;
                load.linearDir = lindir;
                load.linearGbl = lingbl;
                load.linearLcl = linlcl;
                load.tempVar_X = tvx;
                load.tempVar_Y = tvy;
                load.tempVar_Z = tvz;
            end
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class <lelem.html *Lelem*>.
    methods
        %------------------------------------------------------------------
        % Generates element flexural fixed end force (FEF) vector in
        % local xy-plane for an applied linearly distributed load.
        % Output:
        %  fef: 4 position vector with flexural (transversal) FEF's:
        %       fef(1) -> transversal force at initial node
        %       fef(2) -> bending moment at initial node
        %       fef(3) -> transversal force at final node
        %       fef(4) -> bending moment at final node
        function fef = flexuralDistribLoadFEF_XY(load)
            include_constants;
            
            % Initialize transversal load values at end nodes
            qyi = 0;
            qyf = 0;
            
            % Add uniform load contribution in local system
            if ~isempty(load.uniformLcl)
                qyi = qyi + load.uniformLcl(2);
                qyf = qyi;
            end
            
            % Add linear load contribution in local system
            if ~isempty(load.linearLcl)
                qyi = qyi + load.linearLcl(2);
                qyf = qyf + load.linearLcl(5);
            end
            
            % Check if transversal load is not null over element
            if (qyi ~= 0) || (qyf ~= 0)
                L  = load.elem.length;
                L2 = L^2;
                
                % Separate uniform portion from linear portion of tranversal load
                qy0 = qyi;
                qy1 = qyf - qyi;
                
                % Calculate fixed end force vector
                if (load.elem.hingei == CONTINUOUS_END) && (load.elem.hingef == CONTINUOUS_END)
                    
                    fef = [ -(qy0*L/2   + qy1*3*L/20);
                            -(qy0*L2/12 + qy1*L2/30);
                            -(qy0*L/2   + qy1*7*L/20);
                              qy0*L2/12 + qy1*L2/20 ];
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == CONTINUOUS_END)
                    
                    fef = [ -(qy0*3*L/8 + qy1*L/10);
                              0;
                            -(qy0*5*L/8 + qy1*2*L/5);
                              qy0*L2/8  + qy1*L2/15 ];
                    
                elseif (load.elem.hingei == CONTINUOUS_END) && (load.elem.hingef == HINGED_END)
                    
                    fef = [ -(qy0*5*L/8 + qy1*9*L/40);
                            -(qy0*L2/8  + qy1*7*L2/120);
                            -(qy0*3*L/8 + qy1*11*L/40);
                              0 ];
                
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == HINGED_END)
                    
                    fef = [ -(qy0*L/2 + qy1*L/6);
                              0;
                            -(qy0*L/2 + qy1*L/3);
                              0 ];
                end
                
            else
                
                fef = [ 0;
                        0;
                        0;
                        0 ];
            end
        end
        
        %------------------------------------------------------------------
        % Generates element flexural fixed end force (FEF) vector in
        % local xz-plane for an applied linearly distributed load.
        % Output:
        %  fef: 4 position vector with flexural (transversal) FEF's:
        %       fef(1) -> transversal force at initial node
        %       fef(2) -> bending moment at initial node
        %       fef(3) -> transversal force at final node
        %       fef(4) -> bending moment at final node
        function fef = flexuralDistribLoadFEF_XZ(load)
            include_constants;
            
            % Initialize transversal load values at end nodes
            qzi = 0;
            qzf = 0;

            % Add uniform load contribution in local system
            if ~isempty(load.uniformLcl)
                qzi = qzi + load.uniformLcl(3);
                qzf = qzi;
            end
            
            % Add linear load contribution in local system
            if ~isempty(load.linearLcl)
                qzi = qzi + load.linearLcl(3);
                qzf = qzf + load.linearLcl(6);
            end
            
            % Check if transversal load is not null over element
            if (qzi ~= 0) || (qzf ~= 0)
                L  = load.elem.length;
                L2 = L^2;
                
                % Separate uniform portion from linear portion of tranversal load
                qz0 = qzi; 
                qz1 = qzf - qzi;
                
                % Calculate fixed end force vector
                if (load.elem.hingei == CONTINUOUS_END) && (load.elem.hingef == CONTINUOUS_END)
                    
                    fef = [ -(qz0*L/2   + qz1*3*L/20);
                              qz0*L2/12 + qz1*L2/30;
                            -(qz0*L/2   + qz1*7*L/20);
                            -(qz0*L2/12 + qz1*L2/20) ];
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == CONTINUOUS_END)
                    
                    fef = [ -(qz0*3*L/8 + qz1*L/10);
                              0;
                            -(qz0*5*L/8 + qz1*2*L/5);
                            -(qz0*L2/8  + qz1*L2/15) ];
                    
                elseif (load.elem.hingei == CONTINUOUS_END) && (load.elem.hingef == HINGED_END)
                    
                    fef = [ -(qz0*5*L/8 + qz1*9*L/40);
                              qz0*L2/8  + qz1*7*L2/120;
                            -(qz0*3*L/8 + qz1*11*L/40);
                              0 ];
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == HINGED_END)
                    
                    fef = [ -(qz0*L/2 + qz1*L/6);
                              0;
                            -(qz0*L/2 + qz1*L/3);
                              0 ];
                end
                
            else

                fef = [ 0;
                        0;
                        0;
                        0 ];
            end
        end
        
        %------------------------------------------------------------------
        % Generates element flexural fixed end force (FEF) vector in
        % local xy-plane for an applied thermal load.
        % Output:
        %  fef: 4 position vector with flexural (transversal) FEF's:
        %       fef(1) -> transversal force at initial node
        %       fef(2) -> bending moment at initial node
        %       fef(3) -> transversal force at final node
        %       fef(4) -> bending moment at final node
        function fef = flexuralThermalLoadFEF_XY(load)
            include_constants;
            
            % Get temperature gradient relative to element local y-axis
            dty = load.tempVar_Y;
            
            % Check if temperature gradient is not null
            if dty ~= 0
                E     = load.elem.material.elasticity;
                alpha = load.elem.material.thermExp;
                I     = load.elem.section.inertia_z;
                h     = load.elem.section.height_y;
                L     = load.elem.length;
                
                EI = E * I;
                
                % Unitary dimensionless temperature gradient
                tg = (alpha * dty) / h;
                
                % Calculate fixed end force vector
                if (load.elem.hingei == CONTINUOUS_END) && (load.elem.hingef == CONTINUOUS_END)
                    
                    fef = [ 0;
                            tg * (EI);
                            0;
                            tg * (-EI) ];
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == CONTINUOUS_END)
                    
                    fef = [ tg * (-3*EI/(2*L));
                            0;
                            tg * (3*EI/(2*L));
                            tg * (-3*EI/2) ];
                    
                elseif (load.elem.hingei == CONTINUOUS_END) && (load.elem.hingef == HINGED_END)
                    
                    fef = [ tg * (3*EI/(2*L));
                            tg * (3*EI/2);
                            tg * (-3*EI/(2*L));
                            0 ];
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == HINGED_END)
                    
                    fef = [ 0;
                            0;
                            0;
                            0 ];
                end
                
            else
                
                fef = [ 0;
                        0;
                        0;
                        0 ];
            end
        end
        
        %------------------------------------------------------------------
        % Generates element flexural fixed end force (FEF) vector in
        % local xz-plane for an applied thermal load.
        % Output:
        %  fef: 4 position vector with flexural (transversal) FEF's:
        %       fef(1) -> transversal force at initial node
        %       fef(2) -> bending moment at initial node
        %       fef(3) -> transversal force at final node
        %       fef(4) -> bending moment at final node
        function fef = flexuralThermalLoadFEF_XZ(load)
            include_constants;
            
            % Get temperature gradient relative to element local z-axis
            dtz = load.tempVar_Z;
            
            % Check if temperature gradient is not null
            if dtz ~= 0
                E     = load.elem.material.elasticity;
                alpha = load.elem.material.thermExp;
                I     = load.elem.section.inertia_y;
                h     = load.elem.section.height_z;
                L     = load.elem.length;
                
                EI = E * I;
                
                % Unitary dimensionless temperature gradient
                tg = (alpha * dtz) / h;
                
                % Calculate fixed end force vector
                if (load.elem.hingei == CONTINUOUS_END) && (load.elem.hingef == CONTINUOUS_END)
                    
                    fef = [ 0;
                            tg * (-EI);
                            0;
                            tg * (EI) ];
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == CONTINUOUS_END)
                    
                    fef = [ tg * (-3*EI/(2*L));
                            0;
                            tg * (3*EI/(2*L));
                            tg * (3*EI/2) ];
                    
                elseif (load.elem.hingei == CONTINUOUS_END) && (load.elem.hingef == HINGED_END)
                    
                    fef = [ tg * (3*EI/(2*L));
                            tg * (-3*EI/2);
                            tg * (-3*EI/(2*L));
                            0 ];
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == HINGED_END)
                    
                    fef = [ 0;
                            0;
                            0;
                            0 ];
                end
                
            else
                
                fef = [ 0;
                        0;
                        0;
                        0 ];
            end
        end
        
        %------------------------------------------------------------------
        % Computes element transversal displacement in local xy-plane at
        % given cross-section position for an applied linearly distributed load.
        % Output:
        %  dv: transversal displacement in local y-axis
        % Input arguments:
        %  x: cross-section position in element local x-axis
        function dv = flexuralDistribLoadDispl_XY(load,x)
            include_constants;
            
            % Initialize transversal load values at end nodes
            qyi = 0;
            qyf = 0;
            
            % Add uniform load contribution in local system
            if ~isempty(load.uniformLcl)
                qyi = qyi + load.uniformLcl(2);
                qyf = qyi;
            end
            
            % Add linear load contribution in local system
            if ~isempty(load.linearLcl)
                qyi = qyi + load.linearLcl(2);
                qyf = qyf + load.linearLcl(5);
            end
            
            % Check if transversal load is not null over element
            if (qyi ~= 0) || (qyf ~= 0)
                E = load.elem.material.elasticity;
                I = load.elem.section.inertia_z;
                L = load.elem.length;
                
                x2 = x^2;
                x3 = x2 * x;
                x4 = x3 * x;
                x5 = x4 * x;
                L2 = L  * L;
                L3 = L2 * L;
                EI = E  * I;

                % Separate uniform portion from linear portion of tranversal load
                qy0 = qyi;
                qy1 = qyf - qyi;
            
                % Calculate transversal displacements from uniform load portion (vq0)
                % and from linear load portion (vq1)
                if (load.elem.hingei == CONTINUOUS_END) && (load.elem.hingef == CONTINUOUS_END)
                    
                    vq0 = qy0/EI * (L2*x2/24 - L*x3/12 + x4/24);
                    vq1 = qy1/EI * (L2*x2/60 - L*x3/40 + x5/(120*L));
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == CONTINUOUS_END)
                    
                    vq0 = qy0/EI * (L3*x/48  - L*x3/16 + x4/24);
                    vq1 = qy1/EI * (L3*x/120 - L*x3/60 + x5/(120*L));
                    
                elseif (load.elem.hingei == CONTINUOUS_END) && (load.elem.hingef == HINGED_END)
                    
                    vq0 = qy0/EI * (L2*x2/16    - 5*L*x3/48 + x4/24);
                    vq1 = qy1/EI * (7*L2*x2/240 - 3*L*x3/80 + x5/(120*L));
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == HINGED_END)
                    
                    vq0 = qy0/EI * (L3*x/24    - L*x3/12 + x4/24);
                    vq1 = qy1/EI * (7*L3*x/360 - L*x3/36 + x5/(120*L));
                    
                end
                
                dv = vq0 + vq1;
                
            else
                dv = 0;
            end
        end
        
        %------------------------------------------------------------------
        % Computes element transversal displacement in local xz-plane at
        % given cross-section position for an applied linearly distributed load.
        % Output:
        %  dw: transversal displacement in local z-axis
        % Input arguments:
        %  x: cross-section position in element local x-axis
        function dw = flexuralDistribLoadDispl_XZ(load,x)
            include_constants;
            
            % Initialize transversal load value at end nodes
            qzi = 0;
            qzf = 0;
            
            % Add uniform load contribution in local system
            if ~isempty(load.uniformLcl)
                qzi = qzi + load.uniformLcl(3);
                qzf = qzi;
            end
            
            % Add linear load contribution in local system
            if ~isempty(load.linearLcl)
                qzi = qzi + load.linearLcl(3);
                qzf = qzf + load.linearLcl(6);
            end
            
            % Check if transversal load is not null over element
            if (qzi ~= 0) || (qzf ~= 0)
                E = load.elem.material.elasticity;
                I = load.elem.section.inertia_y;
                L = load.elem.length;

                x2 = x^2;
                x3 = x2 * x;
                x4 = x3 * x;
                x5 = x4 * x;
                L2 = L  * L;
                L3 = L2 * L;
                EI = E  * I;

                % Separate uniform portion from linear portion of tranversal load
                qz0 = qzi;
                qz1 = qzf - qzi;
                
                % Calculate transversal displacements from uniform load portion (wq0)
                % and from linear load portion (wq1)
                if (load.elem.hingei == CONTINUOUS_END) && (load.elem.hingef == CONTINUOUS_END)
                    
                    wq0 = qz0/EI * (L2*x2/24 - L*x3/12 + x4/24);
                    wq1 = qz1/EI * (L2*x2/60 - L*x3/40 + x5/(120*L));
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == CONTINUOUS_END)
                    
                    wq0 = qz0/EI * (L3*x/48  - L*x3/16 + x4/24);
                    wq1 = qz1/EI * (L3*x/120 - L*x3/60 + x5/(120*L));
                    
                elseif (load.elem.hingei == CONTINUOUS_END) && (load.elem.hingef == HINGED_END)
                    
                    wq0 = qz0/EI * (L2*x2/16    - 5*L*x3/48 + x4/24);
                    wq1 = qz1/EI * (7*L2*x2/240 - 3*L*x3/80 + x5/(120*L));
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == HINGED_END)
                    
                    wq0 = qz0/EI * (L3*x/24    - L*x3/12 + x4/24);
                    wq1 = qz1/EI * (7*L3*x/360 - L*x3/36 + x5/(120*L));
                    
                end
                
                dw = wq0 + wq1;
                
            else
                dw = 0;
            end
        end
        
        %------------------------------------------------------------------
        % Computes element transversal displacement in local xy-plane at
        % given cross-section position for an applied thermal load.
        % Output:
        %  dv: transversal displacement in local y-axis
        % Input arguments:
        %  x: cross-section position in element local x-axis
        function dv = flexuralThermalLoadDispl_XY(load,x)
            include_constants;
            
            % Get temperature gradient relative to element local y-axis
            dty = load.tempVar_Y;
            
            % Check if temperature gradient is not null
            if dty ~= 0
                alpha = load.elem.material.thermExp;
                h     = load.elem.section.height_y;
                L     = load.elem.length;
                
                x2 = x^2;
                x3 = x2 * x;
                
                % Unitary dimensionless temperature gradient
                tg = (alpha * dty) / h;
                
                % Calculate transversal displacement
                if (load.elem.hingei == CONTINUOUS_END) && (load.elem.hingef == CONTINUOUS_END)
                    
                    dv = 0;
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == CONTINUOUS_END)
                    
                    dv = tg * (-L*x/4 + x2/2 - x3/(4*L));
                    
                elseif (load.elem.hingei == CONTINUOUS_END) && (load.elem.hingef == HINGED_END)
                    
                    dv = tg * (-x2/4 + x3/(4*L));
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == HINGED_END)
                    
                    dv = tg * (-L*x/2 + x2/2);
                    
                end
                
            else
                dv = 0;
            end
        end
        
        %------------------------------------------------------------------
        % Computes element transversal displacement in local xz-plane at
        % given cross-section position for an applied thermal load.
        % Output:
        %  dw: transversal displacement in local z-axis
        % Input arguments:
        %  x: cross-section position in element local x-axis
        function dw = flexuralThermalLoadDispl_XZ(load,x)
            include_constants;
            
            % Get temperature gradient relative to element local z-axis
            dtz = load.tempVar_Z;
            
            % Check if temperature gradient is not null
            if dtz ~= 0
                alpha = load.elem.material.thermExp;
                h     = load.elem.section.height_z;
                L     = load.elem.length;
                
                x2 = x^2;
                x3 = x2 * x;
                
                % Unitary dimensionless temperature gradient
                tg = (alpha * dtz) / h;
                
                % Calculate transversal displacement
                if (load.elem.hingei == CONTINUOUS_END) && (load.elem.hingef == CONTINUOUS_END)
                    
                    dw = 0;
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == CONTINUOUS_END)
                    
                    dw = tg * (-L*x/4 + x2/2 - x3/(4*L));
                    
                elseif (load.elem.hingei == CONTINUOUS_END) && (load.elem.hingef == HINGED_END)
                    
                    dw = tg * (-x2/4 + x3/(4*L));
                    
                elseif (load.elem.hingei == HINGED_END) && (load.elem.hingef == HINGED_END)
                    
                    dw = tg * (-L*x/2 + x2/2);
                    
                end
                
            else
                dw = 0;
            end
        end
    end
end