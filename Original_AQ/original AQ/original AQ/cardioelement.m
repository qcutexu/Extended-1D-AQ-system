% input: iel = element number
%        Asol
%        Qsol = vector of coefficients for the basis function
% 
% output: intA = RHS for A equation without Flux term
%         intQ = RHS for Q equation without Flux term
%
% description: this function returns the an approximate value for
%              the integrals over a particular element.

function [intA, intQ] = cardioelement(iel, Asol, Qsol, variables)

    cardiogrid = variables.cardiogrid;
    pdeg = variables.pdeg;
    xg = variables.xg;
    wg = variables.wg;
    alpha = variables.alpha;
    G0_ = variables.G0_;
    rho = variables.rho;
    type = variables.type;
    nu = variables.nu;
    
    locdim = pdeg + 1;
    
    % obtain values of basis functions at each Gaussian integration point
    [basephi, basephiD, De, xx] = DGelemcalc(pdeg, xg, cardiogrid(iel), cardiogrid(iel+1));

    % get values of Q, A, R0, partialR0, A0 
    [Q, QD] = basiseval(iel, locdim, basephi, basephiD, Qsol);
    [A, AD] = basiseval(iel, locdim, basephi, basephiD, Asol);
    [R0, partialR0, partial2R0, A0] = Stenosis(xx, type, variables);
    
    % calculate values of partialP2 
    %partialP2 = nu*(alpha/(alpha-1))*(Q./A./R0.*partial2R0 + ...
    %            QD./A./R0.*partialR0 - Q.*AD./(A.^2)./R0.*partialR0 - Q./A./A0.*(partialR0.^2));

    % calculate values of F(U) & S(U) at each Gaussian integration point
    FA = Q;
    FQ = alpha*(Q.^2)./A ...
       + G0_/3/rho*(A.^1.5)/(0.18^2);
    %FQ = (alpha- 2/35*(partialR0.^2)).*(Q.^2)./A ...
       %+ G0_/3/rho*(A.^1.5)/(0.18^2);
    

    SA = 0;
    SQ = -2*nu*(alpha/(alpha-1))*(Q./A)...
      + G0_/rho/(0.18^2).*A.*partialR0;
      %ABOVE TERM should not be included in orginal AQ model 
      %- A.*partialP2...
       %+ Q.^2./A.*(-4/35*partialR0.*partial2R0);

    % Gauss quadrature
    intFA = De * sum(wg.*FA.*basephiD, 2);
    intSA = SA;

    intFQ = De * sum(wg.*FQ.*basephiD, 2);
    intSQ = De * sum(wg.*SQ.*basephi, 2);

    intA = intFA + intSA;
    intQ = intFQ + intSQ;

    return