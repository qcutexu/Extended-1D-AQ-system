% input: pdeg = polynomial degree
%        iel = element number
%        which = 1 : A
%              = 2 : Q
% 
% outout: val0 = coefficients for basis polynomials on the interval iel
% 
% description: returns the coefficients for the polynomila basis 
%              for a given initial condition.

function val0 = DGL2proj(pdeg, iel, which, variables)
    
    cardiogrid = variables.cardiogrid;
    wg = variables.wg;
    xg = variables.xg;

    % global cardiogrid wg xg;

    % obtain values of basis functions at each Gaussian integration point
    [basephi, ~, De, xx] = DGelemcalc(pdeg, xg, cardiogrid(iel), cardiogrid(iel + 1));
        
    % obtian initial values at each Gaussian intergration point
    iv = cardioinitial(xx, which, variables);

    % Gauss quadrature
    val0 = De * sum(wg.*iv.*basephi, 2);

    % obtain coefficients for the polynomila basis
    val0 = val0 / De;

    return





