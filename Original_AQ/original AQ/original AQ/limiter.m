% input: Asol
%        Qsol = vector of coefficients for the basis function
%        variables = necessary variables
%
% output: newA
%         newQ = rectified Asol & Qsol using limiter
%
% description: returns rectified Asol & Qsol using limiter

function [newA, newQ] = limiter(Asol, Qsol, variables)
    
    pdeg = variables.pdeg;
    cardiogrid = variables.cardiogrid;
    M = variables.M;
    N = variables.N;
    locdim = pdeg + 1;

    newA = zeros(locdim, N);
    newQ = zeros(locdim, N);


    h = variables.L / N; % grid size
    
    for iel = 1 : N
        if iel == 1 || iel == N
            newA(:, iel) = Asol(:, iel);
            newQ(:, iel) = Qsol(:, iel);
    
        else
            [basephi, basephiD, ~, ~] = DGelemcalc(pdeg, [-1, 1], cardiogrid(iel), cardiogrid(iel+1));
            [A, ~] = basiseval(iel, locdim, basephi, basephiD, Asol);
            [Q, ~] = basiseval(iel, locdim, basephi, basephiD, Qsol);
            
            % average value in iel
            Abar = Asol(1, iel) * basephi(1, 1); 
            Qbar = Qsol(1, iel) * basephi(1, 1);
            
            % average value in iel-1
            Abar1 = Asol(1, iel-1) * basephi(1, 1); 
            Qbar1 = Qsol(1, iel-1) * basephi(1, 1); 
            
            % average value in iel+1
            Abar2 = Asol(1, iel+1) * basephi(1, 1); 
            Qbar2 = Qsol(1, iel+1) * basephi(1, 1); 
    
            DeltaA1 = Abar - Abar1; % DeltaA-
            DeltaQ1 = Qbar - Qbar1; % DeltaQ-
            DeltaA2 = Abar2 - Abar; % DeltaA+
            DeltaQ2 = Qbar2 - Qbar; % DeltaQ+
    
            Ahat = Abar - A(1); % left
            Qhat = Qbar - Q(1);
            Atilde = A(2) - Abar; % right
            Qtilde = Q(2) - Qbar;
    
            Ahat_mod = minmod2(Ahat, DeltaA1, DeltaA2, M, h);
            Atilde_mod = minmod2(Atilde, DeltaA1, DeltaA2, M, h);
            Qhat_mod = minmod2(Qhat, DeltaQ1, DeltaQ2, M, h);
            Qtilde_mod = minmod2(Qtilde, DeltaQ1, DeltaQ2, M, h);
    
            bA = [Asol(1, iel); Abar - Ahat_mod; Abar + Atilde_mod];
            bQ = [Qsol(1, iel); Qbar - Qhat_mod; Qbar + Qtilde_mod];
    
            GA = [1, 0; basephi(:, 1)'; basephi(:, 2)'];
            GQ = [1, 0; basephi(:, 1)'; basephi(:, 2)'];
    
            newA(:, iel) = GA\bA;
            newQ(:, iel) = GQ\bQ;
        end

    end

    return

