cur(n1,n2,val) % adds the stamp of a current source
% (current flowing from n1 to n2).

vol(n1,n2,val) % adds the stamp of a voltage source (n1 positive)
res(n1,n2,val) % adds the stamp of a resistor between node n1 and n2.
ind(n1,n2,val) % adds the inductor stamp.
cap(n1,n2,val) % adds the capacitor stamp.

vccs(nd1,nd2,ni1,ni2,val)
% Stamp of a voltage controlled current source.
% The dependent nodes are nd1 and nd2 (positive current from nd1 to nd2)
% The independent nodes are ni1 and ni2 (positive voltage at ni1).
% Ind1 to Ind2 = val*(Vni1 - Vni2)

vcvs(nd1,nd2,ni1,ni2,val)
% Stamp of a voltage controlled voltage source.
% The dependent nodes are nd1 and nd2 (positive voltage at nd1)
% The independent nodes are ni1 and ni2 (positive voltage at ni1).
% Vnd1 - Vnd2 = val*(Vni1 - Vni2)