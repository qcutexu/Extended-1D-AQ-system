% input:       iel = element number
%              locdim = pdeg + 1
%              basephi = basis function evaluated at a particular point
%              basephiD = derivatives of basis functions evaluated at xloc
%              ysol = vector of coefficients for the basis function
%
% output:      val = value of function
%              valD = value of derivative
%
% description: evaluate function expanding as a linear combo of
%              basis functions at a point xeval

function [val, valD] = basiseval(iel, locdim, basephi, basephiD, ysol)
    
    n = size(basephi, 2);
    ii = (iel - 1)*locdim;
    
    val = zeros(1, n);
    valD = zeros(1, n);
    for jj = 1 : locdim
      val = val + basephi(jj, :)*ysol(ii + jj);
      valD = valD + basephiD(jj)*ysol(ii + jj);
    end
    
    return
