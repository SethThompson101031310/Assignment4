function vccs(nd1,nd2,ni1,ni2,val)
% vccs.m
% Adds stamp for a voltage-controlled current-source to the global circuit representation
% 
%
%   ni1 -------o+          |----------o nd1
%                          |
%                         /|\
%                        / | \    Ind1 to Ind2 = val*(Vni1 - Vni2)
%                        \ | /
%                         \V/ 
%                          |
%   ni2 -------o-          |----------o nd2
%
% The dependent nodes are nd1 and nd2 (positive current from nd1 to nd2)
% The independent nodes are ni1 and ni2 (positive voltage at ni1).
% Ind1 to Ind2 = val*(Vni1 - Vni2)
% Name: Seth Thompson
% Student Number: 101031310

% Defining global variables
global G;

% nd1 and nd2 represent rows, ni1 and ni2 represent columns.
% Check to see which of the four are zeros and place them into the matrix G
% as needed.
if(nd1 ~= 0) && (ni1 ~= 0)
    G(nd1,ni1) = G(nd1,ni1) + val;
end

if(nd1 ~= 0) && (ni2 ~= 0)
    G(nd1,ni2) = G(nd1,ni2) - val;
end

if(nd2 ~= 0) && (ni1 ~= 0)
    G(nd2,ni1) = G(nd2,ni1) - val;
end

if(nd2 ~= 0) && (ni2 ~= 0)
    G(nd2,ni2) = G(nd2,ni2) + val;
end

end

