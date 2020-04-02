function ind(n1,n2,val)
% ind.m
%   Adds stamp for an inductor to the global C-Matrix in circuit representation!
%
% ind(n1,n2,val):
%
%           L = val
%        n1 o---mmmmmm---o n2
% ELEC 4506 Lab 2
% Name: Seth Thompson
% Student Number: 101031310

% Defining global variables.
global G b C;

currentsize = size(G,1); % Obtain the size of the current matrix.
newsize = currentsize + 1; % Index of the new variable to solve.
b(newsize) = 0; % add new row.
G(newsize,newsize) = 0; % add new row/column to G.
C(newsize,newsize) = 0; % add new row/column to C.

% If the first node is not the ground, make the new entries at n1 = 1.
if (n1 ~= 0)
    G(n1,newsize) = 1;
    G(newsize,n1) = 1;
end

% If the second node is not the ground, make the new entries at n2 = -1.
if (n2 ~= 0)
    G(n2,newsize) = -1;
    G(newsize,n2) = -1;
end

% Regardless of what nodes the inductor is connected to, it will still have
% to have its negatie inducance in the diagonal of C
C(newsize,newsize) = -val;

end

