%% ELEC 4700 ASSIGNMENT 4 CODE Q 1-4| Circuit Without Noise

% Name: Seth Thompson
% Student Number: 101031310

close all
clear
clc

% Part 1: Doing a voltage sweep from 0.1V to 10V and obtaining the IV
% characteristics

% To do a voltage sweep of the original system, a for loop must be made to
% simulate the system with n amount of voltage drops across the device. to
% save time, 5 voltages will be tested.

voltageVector = linspace(0.1,10,5);

for SolNumber = 1:numel(voltageVector)

    % Defining the length and width of the box along with Vo

    Length = 2;
    Width = 1;
    Vo = voltageVector(SolNumber);
    % Defining the values for sigma both inside and outside of the boxes.

    sigmaIN = 1e-2;
    sigmaOUT = 1;

    % Defining the dimensions of each of the boxes (wb and lb in figure 3)

    WidthBox = 0.4;
    LengthBox = 0.4;

    % Defining the number of elements that will be in each part of the matrix

    nx = 100*Length;
    ny = 100*Width;

    % Defining the conductivity map

    conductivity = zeros(ny,nx);

    for k = 1:ny
        for l = 1:nx
            % If the element being accesed in the conductivity matrix
            % IS within one of the boxes, set its value to the lower
            % conductivity value
            if(l >= nx*WidthBox && l <= nx-nx*WidthBox && (k >= ny-ny*LengthBox || k <= ny*LengthBox))
                conductivity(k,l) = sigmaIN;
            % Else, put in the higher value
            else
                conductivity(k,l) = sigmaOUT;
            end
        end
    end

    G = sparse(nx*ny,nx*ny);
    B = zeros(nx*ny,1);

    % Populating the G matrix
    for l = 1:nx
        for k = 1:ny

            % Node mapping to put entries into the correct place
            n = k + (l - 1)*ny;

            % Calculating deltas in the x and y direction
            nxm = k + (l - 2)*ny;
            nxp = k + l*ny;
            nym = (k - 1) + (l - 1)*ny;
            nyp = (k + 1) + (l - 1)*ny;

            % Begin inputting all of the correct entires!
            if(l == 1 || l == nx) % Left or right side of the region, set entries to 1
                G(n,n) = 1;
            elseif (k == 1) % We are along the bottom of the region apply finite difference method as needed

                entryYup = (conductivity(k,l)+conductivity(k+1,l))/2;
                entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
                entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;

                G(n,n) = -(entryYup + entryXup + entryXdown);
                G(n,nyp) = entryYup;
                G(n,nxp) = entryXup;
                G(n,nxm) = entryXdown;

            elseif (k == ny) % We are along the top of the region

                entryYdown = (conductivity(k,l)+conductivity(k-1,l))/2;
                entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
                entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;

                G(n,n) = -(entryYdown + entryXup + entryXdown);
                G(n,nym) = entryYdown;
                G(n,nxp) = entryXup;
                G(n,nxm) = entryXdown;
            else % else, apply finite differnce as needed without worrying about going out of bounds...

                % Storing elements from conductivity matrix that will be mapped
                % into the G matrix
                entryYup = (conductivity(k,l)+conductivity(k+1,l))/2;
                entryYdown = (conductivity(k,l)+conductivity(k-1,l))/2;
                entryXup = (conductivity(k,l)+conductivity(k,l+1))/2;
                entryXdown = (conductivity(k,l)+conductivity(k,l-1))/2;

                G(n,n) = -(entryYup + entryYdown + entryXup + entryXdown);
                G(n,nyp) = entryYup;
                G(n,nym) = entryYdown;
                G(n,nxp) = entryXup;
                G(n,nxm) = entryXdown;

            end
        end
    end

    % Populating the B vector next...
    for l = 1:nx
        for k = 1:ny
            % Node mapping to put entries into the correct place
            n = k + (l - 1)*ny;

            % Are we along the left side? if so set the value to Vo
            if (l == 1) 
                B(n) = Vo;
            end

            % Anywhere also it should be zero, but that was defined by using
            % the zeros function to make the vector.
        end
    end

    % Obtaining the solution, V, from V = G\B

    V = G\B;

    % Moing the Solution V from the signle vector to a matrix 

    for l = 1:nx
        for k = 1:ny
            % Node mapping to put entries into the correct place
            n = k + (l - 1)*ny;

            MatrixV(k,l) = V(n);
        end
    end

    % Using the Gradient Function to obtain the electric field from the voltage
    % (E = -del(V))

    [Ex,Ey] = gradient(-MatrixV,1e-9);
    
    % Calculate X and Y component of current density, and then use
    % Pythagoras' Therom to obtain the total current
    currentX(SolNumber) = mean(conductivity(:,1).*Ex(:,1));
    currentY(SolNumber) = mean(conductivity(:,1).*Ex(:,1));
    currentden = sqrt(currentX.^2 + currentY.^2);
       
    % NEW CODE HERE
    
    % the currentden value currently has units of A/m^2, so to get a value
    % that represents the current, the currentden vecotr should be
    % multiplied by the surface area of the face of the component.
    
    % There is one problem with this, there is no height 'z' component
    % given for the device, In this case, the current density will just be
    % multiplied by 100nm, which is the sidelength of the device.
    
    currentVector = currentden.*(100e-9);

end

% With both a voltage vector and a current vector made, plot the data
% points and make a line of best fit!

figure(1)
plot(voltageVector,currentVector,'rx')
xlabel('Voltage (V)')
ylabel('Current (A)')
title({'I-V Characteristics of Device','Seth Thompson | 101031310'})
grid on
hold on
% Creating the line of best fit by making a first order polynomial fit the
% line
coeffs = polyfit(voltageVector,currentVector,1);
yvals = polyval(coeffs,voltageVector);
plot(voltageVector,yvals)
legend('Data Points','Line of Best Fit')
caption = sprintf('y = %fx %f',coeffs(1),coeffs(2));
text(2,2.5,caption)

% From looking at the equation of the line of best fit, it's y-intercept is
% apprxoimetly at 0, and it has a slope of 0.324935A/V. Taking the inverse
% of this value will give the resistance of the device!

fprintf('The value of R3 is %f ohms.\n',1/coeffs(1))

% There will be a seperate file/files reporting on what I did in PA 9

% To simulate the circuit, function files from ELEC 4506 will be sued to
% make the G and C matricies along with the b vector.

global G C b;

G = zeros(5,5); % 5 nodes in the matrix, so initialize G and C as 5x5 
C = zeros(5,5); % zeros matricies and b as a 5x1 solution vector of zeros
b = zeros(5,1);

% This circuit is the same as the one in the MNA PA, with the exception of
% R3's value changing. For this simulation, the new value of R3 will be
% used.

vol(1,0,0.1); % Use 0.1 as a starting point to switch out later.
res(1,2,1);
cap(1,2,0.25);
res(2,0,2);
ind(2,3,0.2);
res(3,0,1/coeffs(1));
vcvs(4,0,3,0,100/(1/coeffs(1)));
res(4,5,0.1);
res(5,0,1000);

% Outputting the G and C matricies
fprintf('The G and C matricies will be output next.\n')
G
C

% ** (G + sC)V = b **

% For the DC case (sweeping from -10V to 10V), ignore the C matrix and only
% focus on the G (as s will always be 0).

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

figure(2)
plot(DCVector,VoDC)
hold on
plot(DCVector,V3DC)
hold off
grid on
title({'Vo and V3 DC Sweep Plot','Seth Thompson | 101031310'})
xlabel('Vin (V)')
ylabel('Calculated Voltage (V)')

% The results here make sense. There should be some form of amplification
% between V3 and Vo, which is seen from the previous plot

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

% Making a plot for the voltage gain in dB

figure(3)
plot(freqVector,abs(Vout))
grid on
title({'Vo AC Sweep Plot','Seth Thompson | 101031310'})
xlabel('frequency (Hz)')
ylabel('Calculated Voltage (V)')

figure(4)
semilogx(freqVector,20*log10(abs(Vout./Vin)))
grid on
title({'Vo/V1 AC Sweep Plot for Gain (dB)','Seth Thompson | 101031310'})
xlabel('frequency (Hz)')
ylabel('Voltage Gain (dB)')

% Looking at the two plots made previously, the results make sense. Both
% plots for the value of Vout and the voltage gain in dB should approach a
% smaller values as the feequency is increased.

% Creating a plot for the randomly perturbed capacitance values

for i = 1:numSteps
    random = pi + 0.05*randn(1,1);
    
   % Cange the C matrix to account for the new C value
   C(1,1) = random;
   C(1,2) = -random;
   C(2,1) = -random;
   C(2,2) = random;
   
   % Solve for the gain again, assume a freq of 1Hz 
   w = 2*pi*1;
   s = 1j*w;
   
    % This time the C matrix must be considered
    Matrix = G + s*C;
    Solution = Matrix\b;

    Vout(i) = Solution(5);
    V1(i) = Solution(1);
    gain2(i) = Solution(5)/Solution(1); 
    absgain2 = abs(gain2);
end

% Setting the cap values back to what they were
cap(1,2,0.25);
% plotting the histogram for gain 
figure(5)
hist(abs(gain2),25)
xlabel('Voltage Gain Range')
ylabel('Number Of Occurences')
title({'Gain Histogram','Seth Thompson | 101031310'})
grid on

% As can be seen in the previous plot, the gain of the system seems to be
% approximetly around 30V/V at 1 Hz.

% PART 4 - Time Domain Analysis.

% A) By looking at the circuit and the previously received results, this
% appears to be some form of amplifier with other components added. 
% Whether its a common source amplifier, common emmitter, or some other
% form is unknown

% B) As the frequency of the input signal increases, the capacitor will act
% more like a short circuit and the inductor will act more like an open
% circut. The inductor will eventually cut off any signal coming into the
% 'amplification' stage of the device, so at higher frequencies there
% should be a smaller output. In short, this should act like a low-pass
% filter.

% C) The trapezoidal rule should be used here to obtain a solution for the
% system in the time domain. It has the equation of 
% x(n+1) ~= x(n) + h/2(x'(n) + x'(n+1)); if we take the time domain
% response of the system to be Gx(t) + Cx'(t) = b(t), then we can
% approximate a solution by using G(xn + xn+1) + C(x'n + x'n+1) = bn +bn+1.

% Substituting the equation for a time domain solution into the trapezoidal
% rule, the equation (G + 2C/h)xn+1 = (2C/h - G)xn + bn + bn+1. Dividing
% the one side by the other will give a result for xn+1 and from there the
% system can be solved!

% D)
% If the circuit is to be simulated for 1second for 1000 steps, then h
% should be 1/1000

h = 1/1000;

% next, the step response will be made that transitions from 0 to 1 at
% 0.03s. The contents of this can be seen in the function file 'step'
timeVector = linspace(0,1,1000);
for i = 1:1000
    stepResponse(i) = Step(timeVector(i));
end

% Doing the time domain analysis using the trapezoidal rule! Note that the
% inital value for Xn is set be be a zero vector since the system starts
% off with no input signal (0V), and as a result the rest of the system
% will reflect that.

X = zeros(8,1);

% Creating a counter to say when a certain part of the 'plot movie' has
% been reached.
counter = 1;
for i = 1:1000
    bn = [0;0;0;0;0;stepResponse(i);0;0];
    if i == numel(timeVector)
        bnplus1 = [0;0;0;0;0;stepResponse(i);0;0];
    end
    if i ~= numel(timeVector)
        bnplus1 = [0;0;0;0;0;stepResponse(i+1);0;0];
    end
    
    nextX = (G + (2*C)/h)\((((2*C)/h) - G)*X + bn + bnplus1);
    VinSTEP(i) = nextX(1);
    VoutSTEP(i) = nextX(5);
    
    % Plotting the graph as a 'movie'
    if counter == 1
        figure(6)
        plot(timeVector(1:i),VinSTEP(1:i),'b')
        grid on
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        title({'Transient Analysis of the System With a Step Response Input','Seth Thompson | 101031310'})
        hold on
        plot(timeVector(1:i),VoutSTEP(1:i),'r')
        xlim([0 1])
        ylim([0 35])
    elseif counter == 1000
        figure(6)
        clf
        plot(timeVector(1:i),VinSTEP(1:i),'b')
        hold on
        plot(timeVector(1:i),VoutSTEP(1:i),'r')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        title({'Transient Analysis of the System With a Step Response Input','Seth Thompson | 101031310'})
        grid on
        xlim([0 1])
        ylim([0 35])
        legend('Input','Output')
        hold off
    else
        figure(6)
        clf
        plot(timeVector(1:i),VinSTEP(1:i),'b')
        hold on
        plot(timeVector(1:i),VoutSTEP(1:i),'r')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        title({'Transient Analysis of the System With a Step Response Input','Seth Thompson | 101031310'})
        xlim([0 1])
        ylim([0 35])
        grid on
    end
    
    X = nextX;
    counter = counter + 1;
    
end

% Plotting fft and fftshift plots

fouriesSTEPin = abs(fft(VinSTEP,2^11));
fouriesSTEPout = abs(fft(VoutSTEP,2^11));

fouriesSTEPinSHIFT = fftshift(fouriesSTEPin);
fouriesSTEPoutSHIFT = fftshift(fouriesSTEPout);

figure(7)
plot(1:2^11,fouriesSTEPin,'b')
hold on
plot(1:2^11,fouriesSTEPout,'r')
hold off
grid on
title({'FFTof Step Response Input','Seth Thompson | 101031310'})
legend('Input','Output')
xlim([0 2048])

figure(8)
plot(linspace(-2^10,2^10,2^11),fouriesSTEPinSHIFT,'b')
hold on
plot(linspace(-2^10,2^10,2^11),fouriesSTEPoutSHIFT,'r')
hold off
title({'FFTshift of Step Response Input','Seth Thompson | 101031310'})
grid on
legend('Input','Output')
xlim([-1028 1028])
% Looking at the results from the time domain simulation above, 

% Next, the simulation will be done with the specified sine wave!
fconst = 1/0.03; %Hz
X = zeros(8,1);
counter = 1;

for i = 1:1000
    bn = [0;0;0;0;0;sin(2*pi*fconst*timeVector(i));0;0];
    if i == numel(timeVector)
        bnplus1 = [0;0;0;0;0;sin(2*pi*fconst*timeVector(i));0;0];
    end
    if i ~= numel(timeVector)
        bnplus1 = [0;0;0;0;0;sin(2*pi*fconst*timeVector(i+1));0;0];
    end
    
    nextX = (G + (2*C)/h)\((((2*C)/h) - G)*X + bn + bnplus1);
    VinSINE(i) = nextX(1);
    VoutSINE(i) = nextX(5);   
    
    % Plotting the graph as a 'movie'
    if counter == 1
        figure(9)
        plot(timeVector(1:i),VinSINE(1:i),'b')
        grid on
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        title({'Transient Analysis of the System With a Sine Wave Input','Seth Thompson | 101031310'})
        hold on
        plot(timeVector(1:i),VoutSINE(1:i),'r')
        xlim([0 1])
        ylim([-3 5])
    elseif counter == 1000
        figure(9)
        clf
        plot(timeVector(1:i),VinSINE(1:i),'b')
        hold on
        plot(timeVector(1:i),VoutSINE(1:i),'r')
        grid on
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        title({'Transient Analysis of the System With a Sine Wave Input','Seth Thompson | 101031310'})
        hold off
        xlim([0 1])
        ylim([-3 5])
        legend('Input','Output')
    else
        figure(9)
        clf
        plot(timeVector(1:i),VinSINE(1:i),'b')
        hold on
        plot(timeVector(1:i),VoutSINE(1:i),'r')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        title({'Transient Analysis of the System With a Sine Wave Input','Seth Thompson | 101031310'})
        xlim([0 1])
        ylim([-3 5])
        grid on
    end
    
    X = nextX;
    counter = counter + 1;
end

% Once again, this plot makes sense. The output waveform is an amplified
% version of the input waveform, with a phase shift applied to it!

% Decreasing the frequency of the input signal increased the voltage at the
% output node while also removing some of the 'transient like response' of
% the system, the last point is most likely due to the fact that the
% capacitor has enough time to fully chaarge and discharge, along with the
% inductor having anough time to disperse its magnetic field.

% Concerning the input frequency, the opposite happens. The small transient
% repsonse stays but the amplitude of the output voltage decreases. Looking
% at the bode plot from before, this makes sense. As the frequency of the
% system increases the voltage gain should also drop, and the opposite is
% true for h=when the frequency decreases. Seeing this in the simulation
% proves that these results are accurate and make sense.

% Plotting fft and fftshift plots

fouriesSINEin = abs(fft(VinSINE,2^11));
fouriesSINEout = abs(fft(VoutSINE,2^11));

fouriesSINEinSHIFT = fftshift(fouriesSINEin);
fouriesSINEoutSHIFT = fftshift(fouriesSINEout);

figure(10)
plot(1:2^11,fouriesSINEin,'b')
hold on
plot(1:2^11,fouriesSINEout,'r')
hold off
title({'FFT of Sine Wave Input','Seth Thompson | 101031310'})
grid on
legend('Input','Output')
xlim([0 2048])

figure(11)
plot(linspace(-2^10,2^10,2^11),fouriesSINEinSHIFT,'b')
hold on
plot(linspace(-2^10,2^10,2^11),fouriesSINEoutSHIFT,'r')
hold off
title({'FFTshift of Sine Wave Input','Seth Thompson | 101031310'})
grid on
legend('Input','Output')
xlim([-1028 1028])

% Next, the gaussian pulse input will be simulated! The vector made below
% will be used as the gaussin pulse input!
gausPulseVector = (1/13.2980)*normpdf(timeVector,0.15,0.03);
X = zeros(8,1);
counter = 1;
for i = 1:1000
    if(timeVector(i) < 0.06)
    bn = [0;0;0;0;0;0;0;0];
        if i == numel(timeVector)
            bnplus1 = [0;0;0;0;0;0;0;0];
        end
        if i ~= numel(timeVector)
            bnplus1 = [0;0;0;0;0;0;0;0];
        end
    end
    
    if(timeVector(i) >= 0.06)
    bn = [0;0;0;0;0;gausPulseVector(i);0;0];
        if i == numel(timeVector)
            bnplus1 = [0;0;0;0;0;gausPulseVector(i);0;0];
        end
        if i ~= numel(timeVector)
            bnplus1 = [0;0;0;0;0;gausPulseVector(i+1);0;0];
        end
    end
    
    
    nextX = (G + (2*C)/h)\((((2*C)/h) - G)*X + bn + bnplus1);
    VinGAUSS(i) = nextX(1);
    VoutGAUSS(i) = nextX(5);   
    
    % Plotting the graph as a 'movie'
    if counter == 1
        figure(12)
        plot(timeVector(1:i),VinGAUSS(1:i),'b')
        grid on
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        title({'Transient Analysis of the System With a Gaussian Pulse Input','Seth Thompson | 101031310'})
        hold on
        xlim([0 1])
        ylim([-1 20])
        plot(timeVector(1:i),VoutGAUSS(1:i),'r')
    elseif counter == 1000
        figure(12)
        clf
        plot(timeVector(1:i),VinGAUSS(1:i),'b')
        hold on
        plot(timeVector(1:i),VoutGAUSS(1:i),'r')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        title({'Transient Analysis of the System With a Gaussian Pulse Input','Seth Thompson | 101031310'})
        grid on
        xlim([0 1])
        ylim([-1 20])
        legend('Input','Output')
        hold off
    else
        figure(12)
        clf
        plot(timeVector(1:i),VinGAUSS(1:i),'b')
        hold on
        plot(timeVector(1:i),VoutGAUSS(1:i),'r')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        title({'Transient Analysis of the System With a Gaussian Pulse Input','Seth Thompson | 101031310'})
        xlim([0 1])
        ylim([-1 20])
        grid on
    end
    
    X = nextX;
    counter = counter + 1;
end

fouriesGAUSSin = abs(fft(VinGAUSS,2^11));
fouriesGAUSSout = abs(fft(VoutGAUSS,2^11));

fouriesGAUSSinSHIFT = fftshift(fouriesGAUSSin);
fouriesGAUSSoutSHIFT = fftshift(fouriesGAUSSout);

figure(13)
plot(1:2^11,fouriesGAUSSin,'b')
hold on
plot(1:2^11,fouriesGAUSSout,'r')
hold off
title({'FFT of Gaussian Pulse Input','Seth Thompson | 101031310'})
grid on
legend('Input','Output')
xlim([0 2048])

figure(14)
plot(linspace(-2^10,2^10,2^11),fouriesGAUSSinSHIFT,'b')
hold on
plot(linspace(-2^10,2^10,2^11),fouriesGAUSSoutSHIFT,'r')
hold off
title({'FFTshift of Gaussian Pulse Input','Seth Thompson | 101031310'})
grid on
legend('Input','Output')
xlim([-1028 1028])
% With that, the simulation for the gaussian pulse is done!

% Increaseing the time step of each simulation (h) will lead to more
% innacurate results. The trapezoidal rule works because very small steps
% are taken between data points, and in between that step the next data
% point is calculated with a given amount of error. In this case, larger
% values of h tend to give a 'sawtooth like' result, as shown next.

h = 1/100;

% Use the sine wave plot code to show why h should not be too large
X = zeros(8,1);
counter = 1;
timeVectorDummy = linspace(0,1,100);
fconst = 1/0.03; %Hz
for i = 1:100
    bn = [0;0;0;0;0;sin(2*pi*fconst*timeVectorDummy(i));0;0];
    if i == numel(timeVectorDummy)
        bnplus1 = [0;0;0;0;0;sin(2*pi*fconst*timeVectorDummy(i));0;0];
    end
    if i ~= numel(timeVectorDummy)
        bnplus1 = [0;0;0;0;0;sin(2*pi*fconst*timeVectorDummy(i+1));0;0];
    end
    
    nextX = (G + (2*C)/h)\((((2*C)/h) - G)*X + bn + bnplus1);
    VinSINEDummy(i) = nextX(1);
    VoutSINEDummy(i) = nextX(5);   
    
    % Plotting the graph as a 'movie'
    if counter == 1
        figure(15)
        plot(timeVectorDummy(1:i),VinSINEDummy(1:i),'b')
        grid on
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        title({'Transient Analysis With Too Large a Time Step','Seth Thompson | 101031310'})
        hold on
        plot(timeVectorDummy(1:i),VoutSINEDummy(1:i),'r')
        xlim([0 1])
    elseif counter == 1000
        figure(15)
        clf
        plot(timeVectorDummy(1:i),VinSINEDummy(1:i),'b')
        hold on
        plot(timeVectorDummy(1:i),VoutSINEDummy(1:i),'r')
        grid on
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        title({'Transient Analysis With Too Large a Time Step','Seth Thompson | 101031310'})
        hold off
        xlim([0 1])
        legend('Input','Output')
    else
        figure(15)
        clf
        plot(timeVectorDummy(1:i),VinSINEDummy(1:i),'b')
        hold on
        plot(timeVectorDummy(1:i),VoutSINEDummy(1:i),'r')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        title({'Transient Analysis With Too Large a Time Step','Seth Thompson | 101031310'})
        xlim([0 1])
        grid on
    end
    
    X = nextX;
    counter = counter + 1;
end