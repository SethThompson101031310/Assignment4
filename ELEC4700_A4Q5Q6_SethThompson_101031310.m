%% ELEC 4700 ASSIGNMENT 4 CODE Q5 and Q6 | Circuit With Noise

% Name: Seth Thompson
% Student Number: 101031310

close all
clear
clc

% Repeating the analysis to obtain a value for R3...

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

coeffs = polyfit(voltageVector,currentVector,1);
R3 = 1/coeffs(1);

% Creating new G, C, and b variables for the system with noise

global G C b;

G = zeros(5,5); % 5 nodes in the matrix, so initialize G and C as 5x5 
C = zeros(5,5); % zeros matricies and b as a 5x1 solution vector of zeros
b = zeros(5,1);

vol(1,0,0.1); % Use 0.1 as a starting point to switch out later.
res(1,2,1);
cap(1,2,0.25);
res(2,0,2);
ind(2,3,0.2);
res(3,0,1/coeffs(1));
vcvs(4,0,3,0,100/(1/coeffs(1)));
res(4,5,0.1);
res(5,0,1000);
% Adding in a dummy variabel for the new current source.
cur(3,0,0.001);
% Adding in the new capacitor
cap(3,0,0.00001);

% outputting the new C matrix
fprintf('The new C matrix is shown next.\n')
C

numSteps = 2000;
% Set Vin value to 1V
b(6,1) = 1;

% Doing a new time domain analysis with the new capacitor and noise-current
% source.

h = 1/2000;
timeVector = linspace(0,1,2000);
fconst = 1/0.03; %Hz
X = zeros(8,1);

for i = 1:2000
    % Storing the random noise in the randVal variable.
    randVal = -0.001 +(randn)/10000;
    
    bn = [0;0;randVal;0;0;sin(2*pi*fconst*timeVector(i));0;0];
    if i == numel(timeVector)
        bnplus1 = [0;0;randVal;0;0;sin(2*pi*fconst*timeVector(i));0;0];
    end
    if i ~= numel(timeVector)
        bnplus1 = [0;0;randVal;0;0;sin(2*pi*fconst*timeVector(i+1));0;0];
    end
    
    nextX = (G + (2*C)/h)\((((2*C)/h) - G)*X + bn + bnplus1);
    VinSINE(i) = nextX(1);
    VoutSINE(i) = nextX(5);   
    
    % Note that this time, it wasnt asked to show this plotting over time
    
    X = nextX;
end

figure(1)
plot(timeVector,VinSINE,'b')
grid on
xlabel('Time (s)')
ylabel('Voltage (V)')
title({'Transient Analysis of the System With Noise (Cn = 0.00001, h = 1/2000)','Seth Thompson | 101031310'})
hold on
plot(timeVector,VoutSINE,'r')
legend('Input','Output')
xlim([0 1])

% plotting the fft of the output and input
fftIN = abs(fft(VinSINE,2^11));
fftOUT = abs(fft(VoutSINE,2^11));
fftINSHIFT = fftshift(fftIN);
fftOUTSHIFT = fftshift(fftOUT);

figure(2)
plot(1:1:2^11,fftIN,'b')
hold on
plot(1:1:2^11,fftOUT,'r')
hold off
grid on
title({'FFT of Sine Wave Input and Output (Noisy)','Seth Thompson | 101031310'})
legend('Input','Output')
xlim([0 2050])

figure(3)
plot(-2^10:1:2^10-1,fftINSHIFT,'b')
hold on
plot(-2^10:1:2^10-1,fftOUTSHIFT,'r')
hold off
grid on
title({'FFTshift of Sine Wave Input and Output (Noisy)','Seth Thompson | 101031310'})
legend('Input','Output')
xlim([-1250 1250])

% The fourier transform here makes sense, the major peaks are where the
% frequency of the sine wave is!

% Doing the simulation again, but this time with a larger capacitor.
X = zeros(8,1);
C(3,3) = 0.00001*1000;

for i = 1:2000
    % Storing the random noise in the randVal variable.
    randVal = -0.001 +(randn)/10000;
    
    bn = [0;0;randVal;0;0;sin(2*pi*fconst*timeVector(i));0;0];
    if i == numel(timeVector)
        bnplus1 = [0;0;randVal;0;0;sin(2*pi*fconst*timeVector(i));0;0];
    end
    if i ~= numel(timeVector)
        bnplus1 = [0;0;randVal;0;0;sin(2*pi*fconst*timeVector(i+1));0;0];
    end
    
    nextX = (G + (2*C)/h)\((((2*C)/h) - G)*X + bn + bnplus1);
    VinSINE2(i) = nextX(1);
    VoutSINE2(i) = nextX(5);   
    
    % Note that this time, it wasnt asked to show this plotting over time
    
    X = nextX;
end

figure(4)
plot(timeVector,VinSINE2,'b')
grid on
xlabel('Time (s)')
ylabel('Voltage (V)')
title({'Transient Analysis of the System With Noise (Cn = 0.01, h = 1/2000)','Seth Thompson | 101031310'})
hold on
plot(timeVector,VoutSINE2,'r')
legend('Input','Output')
xlim([0 1])

% Doing the simulation again, but this time with a smaller capacitance


% Doing the simulation again, but this time with a larger capacitor.
X = zeros(8,1);
C(3,3) = 0.00001/1000;

for i = 1:2000
    % Storing the random noise in the randVal variable.
    randVal = -0.001 +(randn)/10000;
    
    bn = [0;0;randVal;0;0;sin(2*pi*fconst*timeVector(i));0;0];
    if i == numel(timeVector)
        bnplus1 = [0;0;randVal;0;0;sin(2*pi*fconst*timeVector(i));0;0];
    end
    if i ~= numel(timeVector)
        bnplus1 = [0;0;randVal;0;0;sin(2*pi*fconst*timeVector(i+1));0;0];
    end
    
    nextX = (G + (2*C)/h)\((((2*C)/h) - G)*X + bn + bnplus1);
    VinSINE3(i) = nextX(1);
    VoutSINE3(i) = nextX(5);   
    
    % Note that this time, it wasnt asked to show this plotting over time
    
    X = nextX;
end

figure(5)
plot(timeVector,VinSINE3,'b')
grid on
xlabel('Time (s)')
ylabel('Voltage (V)')
title({'Transient Analysis of the System With Noise (Cn = 1e-8, h = 1/2000)','Seth Thompson | 101031310'})
hold on
plot(timeVector,VoutSINE3,'r')
legend('Input','Output')
xlim([0 1])

% After looking at the three plots with different values of Cn, it becomes
% clear that and the value of Cn approaches zero, the bandwidth of the
% system increases an lets more of the signal through at thigher
% frequencies, and as the value increases, the band width decreases and
% less of the signal is able to go through at higher frequencies.

% Looking at the equation for the impedance of a capacitor, (-j/sC), it
% becomes apperent that as the value of C decreasses, the impedance should
% get large, and as it increases the impedance gets smaller. IF the
% impedance gets smaller, then the capacitor will act more like a short
% circuit and in this case look more like a short to ground, hence letting
% less of a signal control what the dependant source outputs; the opposite
% is true for if the impedance increases! less of the signal will be lost
% and hence more of the signal can control what the dependant source
% outputs! that is why the outputs appear as they do!

% Doing the simulation again, but this larger step size
h = 1/100;
X = zeros(8,1);
C(3,3) = 0.00001;
timeVector = linspace(0,1,100);
for i = 1:100
    % Storing the random noise in the randVal variable.
    randVal = -0.001 +(randn)/10000;
    
    bn = [0;0;randVal;0;0;sin(2*pi*fconst*timeVector(i));0;0];
    if i == numel(timeVector)
        bnplus1 = [0;0;randVal;0;0;sin(2*pi*fconst*timeVector(i));0;0];
    end
    if i ~= numel(timeVector)
        bnplus1 = [0;0;randVal;0;0;sin(2*pi*fconst*timeVector(i+1));0;0];
    end
    
    nextX = (G + (2*C)/h)\((((2*C)/h) - G)*X + bn + bnplus1);
    VinSINE4(i) = nextX(1);
    VoutSINE4(i) = nextX(5);   
    
    % Note that this time, it wasnt asked to show this plotting over time
    
    X = nextX;
end

figure(6)
plot(timeVector,VinSINE4,'b')
grid on
xlabel('Time (s)')
ylabel('Voltage (V)')
title({'Transient Analysis of the System With Noise (Cn = 0.00001, h = 1/100)','Seth Thompson | 101031310'})
hold on
plot(timeVector,VoutSINE4,'r')
legend('Input','Output')
xlim([0 1])

% The result obtained here makes sense too! As mentioned before, having a
% smaller step size will decrease the accuracy of the simulation and give
% less-accurate results. This can be seen by comparing figure 1 to figure
% 4, which should have the same resluts if the same step size was used.

% QUESTION 6 Answer

% A) The new dependant source is now non-linear in fashion, so to simulate
% the system, a newton-rapshon method would need to be used to approximate
% the solution of the system at each time step due to the non-linear
% component.

% B) In general, the new solution to the equation whould have the form
% GV + CV' + F(x) = b, where F(x) is a vector containing all of the
% non-linear elements of the system. It should also be known that the
% equation to the newton-rapshon serach has the form 
% xn+1 = xn - (phi(x)/phi'(x)). In this case, the x(n)'s will be the
% solutions to the voltage vector, and the general soltuion to the systems
% turns into (G + 2C/h)xn+1 + F(xn+1) = (2C/h -G)xn -F(xn) +bn + bn+1. Also
% note that for NR search, a jacobian matrix of the system must be taken to act
% as phi'(x).

% The system can be solved by completing a time step, but then after that
% using the newton rapshon method to converge to the real result. This can
% be done by seeing when the difference between two different iterations of
% the search goes below a threshold value. After that, the next time step
% should be completed, then the NR method again, until the simulation is
% done!.