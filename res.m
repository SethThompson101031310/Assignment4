function res(n1,n2,val)
% res.m
% Adds stamp for a resistor to the global G-Matrix in circuit representation!
%
% res(n1,n2,val):
%
%           R = val
%        n1 o---/\/\/\---o n2
% ELEC 4506 Lab 2
% Name: Seth Thompson
% Student Number: 101031310

% Defining global variables.
global G;

% Is node 1 ground? if not then add the inverse of the resistance. 
if (n1 ~= 0)
    G(n1,n1) = G(n1,n1) + (1/val);
end
% Is node 2 ground? if not then add the inverse of the resistance.
if (n2 ~= 0)
    G(n2,n2) = G(n2,n2) + (1/val);
end

% Are either of the nodes ground? if not then subtract the inverse of the
% resistance from the off-diagonal entries.
if (n1 ~= 0) && (n2 ~= 0)
    G(n1,n2) = G(n1,n2) - (1/val);
    G(n2,n1) = G(n2,n1) - (1/val);
end
end

