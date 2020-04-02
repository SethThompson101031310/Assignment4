function cur(n1,n2,val)
% cur.m:
% Add stamp for current source to the global circuit representation
%            
%               Csrc=val (A)
%       n1 o-----(--->)-----o n2
%               +      -
%                 VolCS
%         
% ELEC4506, Lab 2
% Name: Seth Thompson
% Student Number: 101031310

% Defining Global Varibles
global b;

% Depending on the direction of the current, the b vector will have the
% value of the current source added or subtracted to it.
if(n1 ~= 0)
    b(n1) = b(n1) - val;
end

if(n2 ~= 0)
    b(n2) = b(n2) + val;
end

end

