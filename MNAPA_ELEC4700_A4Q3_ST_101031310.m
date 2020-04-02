% Name: Seth Thompson
% Student Number: 101031310

% ELEC 4700 MNA PA | February 28th

close all
clear
clc

% The differenctial equations of the system from KCL are as follows...

% 0 = (V2-V1)/R1 + (V2' - V1')C + IE
% 0 = (V1-V2)/R1 + (V1' - V2')C + (0-V2)/R2 + (V3 - V2)dt/L
% 0 = (V2-V3)dt/L + (0-V3)/R3
% 0 = IDP + (V5-V4)/R4 
% 0 = (V4-V5)/R4 + (0-V5)Ro

% These equations look like the following in the frequency domain
% 0 = (V2-V1)/R1 + (V2 - V1)sC + IE
% 0 = (V1-V2)/R1 + (V1 - V2)sC + (0-V2)/R2 + (V3 - V2)/sL
% 0 = (V2-V3)/sL + (0-V3)/R3
% 0 = IDP + (V5-V4)/R4 
% 0 = (V4-V5)/R4 + (0-V5)Ro

% Defining global G, C and B

global G C b;

% Starting the system off with matricies the same size as the number of
% nodes, functions from ELEC 4506 will populate each matrix/vector.

G = zeros(5,5);
C = zeros(5,5);
b = zeros(5,1);

% Filling in the matricies/vector

vol(1,0,5) % Use 5 as a dummy variable
res(1,2,1)
cap(1,2,0.25)
res(2,0,2)
ind(2,3,0.2)
res(3,0,10)
res(4,5,0.1)
res(5,0,1000)
% Jury-rig CCVS by dividing gain by R3
vcvs(4,0,3,0,100/10)

% The DC sweep can now be done, must sweep Vin from -10V to 10V

numSteps=200;
DCVector = linspace(-10,10,numSteps);

% Use a loop to obtain all solutions
for i = 1:numSteps
    
    % DC solution, so no need to use C matrix (s = 0)
    b(6,1) = DCVector(i);
    
    % Make a vector that holds the full solution
    Solution = G\b;
    
    % Hold Vo and V3 answers in their own vectors
    VoDC(i) = Solution(5);
    V3DC(i) = Solution(3);
end

% Plotting both Vo and V3 on the same plot

figure(1)
plot(DCVector,VoDC)
hold on
plot(DCVector,V3DC)
hold off
grid on
title({'Vo and V3 DC Sweep Plot','Seth Thompson | 101031310'})
xlabel('Vin (V)')
ylabel('Calculated Voltage (V)')

% Making another loop for the voltage gain

freqVector = linspace(1,100,numSteps);
% Set Vin value to 1V
b(6,1) = 1;
for i = 1:numSteps
   
    % Making the s variable
    w = 2*pi*freqVector(i);
    s = 1j*w;
    
    % This time the C matrix must be considered
    Matrix = G + s*C;
    Solution = Matrix\b;
    
    Vout(i) = Solution(5);
    Vin(i) = Solution(1);
    gain(i) = Solution(5)/Solution(1); 
end

% Making a plot for the output voltage

figure(2)
semilogx(freqVector,abs(Vout))
grid on
title({'Vo AC Sweep Plot','Seth Thompson | 101031310'})
xlabel('frequency (Hz)')
ylabel('Calculated Voltage (V)')

% Creating a plot for the voltage gain in dB
figure(3)
semilogx(freqVector,20*log10(abs(gain)))
grid on
title({'Voltage Gain AC Sweep Plot','Seth Thompson | 101031310'})
xlabel('frequency (Hz)')
ylabel('Voltage Gain (dB)')

% Creating a plot for the randomly perturbed capacitance values

for i = 1:numSteps
    random = pi + 0.05*randn(1,1);
    
   % Cange the C matrix to account for the new C value
   C(1,1) = random;
   C(1,2) = -random;
   C(2,1) = -random;
   C(2,2) = random;
   
   % Store C value in vector
   CStore(i) = random;
   % Solve for the gain again 
   w = 2*pi*freqVector(i);
   s = 1j*w;
   
    % This time the C matrix must be considered
    Matrix = G + s*C;
    Solution = Matrix\b;

    Vout(i) = Solution(5);
    Vin(i) = Solution(1);
    gain(i) = Solution(5)/Solution(1); 
end

% plotting the histogram for gain 
figure(4)
hist(abs(gain),20)
grid on
title('Gain Histogram')