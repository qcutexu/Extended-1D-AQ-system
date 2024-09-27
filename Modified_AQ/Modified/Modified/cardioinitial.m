% input: xx = global point
%        which = 1 : A
%              = 2 : Q
% 
% output: val = value of initial mass flux or area at the point xx
%
% description: returns value of initial condition at the point xx

function val = cardioinitial(xx, which, variables)
    
    type = variables.type;
    % global type
    
    if which == 1 % A
    
        [~, ~, ~, val] = Stenosis(xx, type, variables);
    
    elseif which == 2 % Q
    
        val = 0;
    
    end
    
    return