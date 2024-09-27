% input: pdeg = polynomial degree
%        xloc = local point in reference element [-1,1]
%        xe = global left boundary of element e
%        xe1 = global right boundary of element e
%
% output: basephi = local basis functions (Legendre polynomials)
%         basephiD = derivatives of basis functions evaluated at xloc
%         De = determinant of Jacobian of transformation
%         xx = global coordinate of xloc      
%
% description: returns necessary information for a particular element

function [basephi, basephiD, De, xx] = DGelemcalc(pdeg, xloc, xe, xe1)
    
    n = length(xloc);

    % initialize size of basephi
    basephi = zeros(pdeg + 1, n);
    basephiD = zeros(pdeg + 1, 1);

    % compute global point
    xx = 0.5*(xe1 - xe)*xloc + 0.5*(xe + xe1);
    
    % compute determinant of jacobian
    De = 0.5*(xe1 - xe);
    
    % basis
    basephi(1, :) = 1/sqrt(2); 
    basephi(2, :) = sqrt(1.5) * xloc;
    
    % derivative of basis
    basephiD(1) = 0; 
    basephiD(2) = sqrt(1.5)*(1/De);
    
    return