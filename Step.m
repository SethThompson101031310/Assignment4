function [Voltage] = Step(time)
% Step | Returns the value of the specified step function at a certain time 
%   time: time of the system
%   Voltage: Voltage at the point in time

% Name: Seth Thompson
% Student Number: 101031310

if (time >= 0 && time < 0.03)
    Voltage = 0;
elseif (time >= 0.03)
    Voltage = 1;
end


end

