% input: Asol
%        Qsol = vector of coefficients for the basis function
%
% 
% output: Ain
%         Qin = Inlet boundary values determined by first extrapolating the 
%                left-to-right moving Riemann invariant to the boundary
%         Aout
%         Qout = Outlet boundary values determined by first extrapolating the 
%                right-to-left moving Riemann invariant to the boundary
% 
% description: specify boundary data on boundary nodes (Velocity-driven)

function [Ain, Qin, Aout, Qout] = cardiobdryext_V(Asol, Qsol, variables)
    
    cardiogrid = variables.cardiogrid;
    rho = variables.rho;
    N = variables.N;
    type = variables.type;
    G0_ = variables.G0_;
    pdeg = variables.pdeg;
    dt = variables.dt;

    %global cardiogrid rho N type G0_ pdeg dt

    locdim = pdeg + 1;

    %%%%%% INLET %%%%%%
    iel = 1;
    xloc = -1;
    [basephi, basephiD, ~, xe] = DGelemcalc(pdeg, xloc, cardiogrid(iel), cardiogrid(iel+1));

    % compute A0 at xe
    [~, ~, ~, A0] = Stenosis(xe,type, variables);

    % compute AR, QR at xe
    [AR, ~] = basiseval(iel, locdim, basephi, basephiD, Asol);
    [QR, ~] = basiseval(iel, locdim, basephi, basephiD, Qsol);

    % compute eigenvalues at inlet
    lam2 = (QR/AR) - sqrt(G0_/rho/2*sqrt(AR)/(0.18^2));
    ind = 1;
    while cardiogrid(ind+1) < -dt*lam2
        ind = ind + 1;
    end

    % compute Q and A evaluated at cardiogrid(1) - lam2*dt
    xhatloc = (2/(cardiogrid(ind+1) - cardiogrid(ind)))*(-lam2*dt - 0.5*(cardiogrid(ind+1) + cardiogrid(ind)));
    [basephihat, basephihatD, ~, ~] = DGelemcalc(pdeg, xhatloc, cardiogrid(ind), cardiogrid(ind+1));
    [Ahat, ~] = basiseval(ind, locdim, basephihat, basephihatD, Asol);
    [Qhat, ~] = basiseval(ind, locdim, basephihat, basephihatD, Qsol);

    % compute Riemann Invariant at xhatloc
    w2 = (Qhat/Ahat) - 4*sqrt(G0_/2/rho/(0.18^2))*(Ahat^0.25);

    % specify boundary data on boundary nodes
    Qin = A0 * 22.5;
    fun = @(x) Qin/x - w2 - 4*sqrt(G0_/2/rho/(0.18^2)) * x^0.25;
    Ain = fzero(fun, Ahat);
    % Ain = A0;

    %%%%%% OUTLET %%%%%%
    iel = N;
    xloc = 1;
    [basephi, basephiD, ~, xe1] = DGelemcalc(pdeg, xloc, cardiogrid(iel), cardiogrid(iel+1));

    % compute A0 at xe1
    [~, ~, ~, A0] = Stenosis(xe1,type, variables);

    % compute AL, QL at xe1
    [AL, ~] = basiseval(iel, locdim, basephi, basephiD, Asol);
    [QL, ~] = basiseval(iel, locdim, basephi, basephiD, Qsol);

    % compute eigenvalues at outlet
    lam1 = (QL/AL) + sqrt(G0_/rho/2*sqrt(AL)/(0.18^2));
    ind = N;
    while cardiogrid(ind) > (cardiogrid(N+1) -dt*lam1)
        ind = ind - 1;
    end

    % compute Q and A evaluated at cardiogrid(N+1) - lam1*dt
    xhatloc = (2/(cardiogrid(ind+1) - cardiogrid(ind)))*(cardiogrid(N+1)-lam1*dt - 0.5*(cardiogrid(ind+1) + cardiogrid(ind)));
    [basephihat, basephihatD, ~, ~] = DGelemcalc(pdeg, xhatloc, cardiogrid(ind), cardiogrid(ind+1));
    [Ahat, ~] = basiseval(ind, locdim, basephihat, basephihatD, Asol);
    [Qhat, ~] = basiseval(ind, locdim, basephihat, basephihatD, Qsol);

    % compute Riemann Invariant at xhatloc
    w1 = (Qhat/Ahat) + 4*sqrt(G0_/2/rho/(0.18^2))*Ahat^0.25;
    
    % specify boundary data on boundary nodes
    Aout = A0;
    Qout = Aout*w1 - 4*sqrt(G0_/2/rho/(0.18^2))*Aout^1.25;
    
    return