% input: a = u_tilde
%        b = delta_plus
%        c = delta_minus
%
% output: m = minmod value for (a, b, c)
%
% description: minmod function
% m(a, b ,c) = s*min(|a|, |b|, |c|}, s = sgn(a) = sgn(b) = sgn(c)
%            = 0, otherwise

function m = minmod(a, b, c)
    X = [a, b, c];
    s = sum(sign(X)) / length(X);

    if abs(s) == 1
        m = min(abs(X));
    else
        m = 0;
    end
    
    return
    