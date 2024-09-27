% input: a = u_tilde
%        b = delta_plus
%        c = delta_minus
%        M = coefficient
%        h = grid size
%
% output: m = rectified minmod value for (a, b, c)
%
% description: rectified minmod function
% m(a, b ,c) = a, if |a1| <= Mh^2
%            = minmod(a, b, c), otherwise

function m = minmod2(a, b, c, M, h)
    
    if abs(a) <= M*h^2
        m = a;
    else
        m = minmod(a, b, c);
    end

    return
