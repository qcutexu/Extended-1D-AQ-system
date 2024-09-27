% input: iel = element number
%        Asol
%        Qsol = vector of coefficients for the basis function
%        Bext = specified boundary data
% 
% output: Atmp
%         Qtmp = intermediate stage in RK scheme
% 
% description: intermediate stage in RK scheme

function [Atmp, Qtmp] = RKrhs(iel, Asol, Qsol, Bext, variables)
    
    pdeg = variables.pdeg;
    % obtain Lax-Friedrichs flux values for element iel 
    [fluxA, fluxQ] = LFflux(iel, pdeg, Asol, Qsol, Bext, variables);
    
    % inegral terms over element iel
    [intA, intQ] = cardioelement(iel, Asol, Qsol, variables);

    % compute flux term
    totfluxA = fluxA(:, 2) - fluxA(:, 1);
    totfluxQ = fluxQ(:, 2) - fluxQ(:, 1);
    
    Atmp = -totfluxA + intA;
    Qtmp = -totfluxQ + intQ;

    return