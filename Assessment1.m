%Assessment code written by Kieran Currie-Cathey 9/12/21
clc,clear;




%list of inputs / variables to be used in given equations
L=input('What is the total length of the beam (m)? ');
Lab=input('What is the distance between the two supports (m)? ');
W=input('What is the weight per unit length of the beam (N/m)? ');
num_of_loads=input('What is the total number of loads acting on the beam? ');
D=input('What is the flexural rigidity (Nm^2)? ');
w = zeros(1, length(num_of_loads));
p = zeros(1, length(num_of_loads));




%using a for loop to allow for any number of loads to used on this code,
%changes the number of inputs here based on the integer given for
%num_of_loads, the while loop
%gives an error if the user enters a value for the distance of the given
%load in excess of the distance of the second support, and will continue to ask for an input until the condition is satisifed, then will move onto the next load
for a=1:num_of_loads
    w(a)=input(['Magnitude of load ' int2str(a) ', with posotive being the load applied downwards, (N): ']);
    p(a)=input(['Input distance from left hand support of load ' int2str(a) ' in (m): ']);
    while p(a)>Lab
        disp('ERROR: the distance of the load given must be less than that of the second support, please re-enter the value. ');
        p(a)=input(['Input distance from left hand support of load ' int2str(a) ' in (m): ']);
    end
end

%dividing the beam up into 1000 increments (arbitrary number, the larger
%the more accurate the calculations)
x=linspace(0,L,1000);

%finding all the points of x<=Lab and Lab<x<=L
f1=find(x<=Lab);
f2=find((x>Lab)&(x<=L));

%assinging all the points found by 'find' into vectors for use in equations
x1=x(f1);
x2=x(f2);

%calculating the delfection due to the weight of the beam between the two
%given bounds
Y01=(W/(24*D)).*(-x1.^4+2*(2*L-(L^2/Lab)).*x1.^3+(2*L^2*Lab-4*L*Lab^2+Lab^3).*x1);
Y02=(W/(24*D)).*(-3*Lab^3+8*L*Lab^2-4*L^2*Lab).*(x2-Lab);
Y0=[Y01 Y02];

%dummy variable made to be later used in the for loop / summation
b=0;

%now need to find all the requisite points within the given bounds for each
%equation, and each load, for loop required so calculations can be done for
%each load
%qk must be defined within the for loop as it is dependent on the load in
%question
for k=1:num_of_loads
    point=p(k);
    qk=Lab-point;%defining qk for each iteration
    f3=find(x<=point);%finding all the points of x that satisfy the given bounds in these 3 lines (f3 f4 f5)
    f4=find((x>point & x<=Lab));
    f5=find((x>Lab)&(x<=L));
    x3=x(f3);%creating seperate vectors of all the values found by the find functions so they can be used in vector products for the calculations
    x4=x(f4);
    x5=x(f5);
    Yk1=((w(k)*qk)/(6*D*Lab)).*(x3.^3-(Lab^2-qk^2).*x3);%equations to calculate the deflection in given intervals for each load
    Yk2=((w(k)*qk)/(6*D*Lab))*(x4.^3-(Lab/qk).*(x4-point).^3-(Lab^2-qk^2).*x4);
    Yk3=((w(k)*qk)/(6*D*Lab))*(3*Lab^2-((3*Lab)/qk).*(Lab-point)^2-(Lab^2-qk^2)).*(x5-Lab);
    Yk=[Yk1 Yk2 Yk3];%creating a vector containing all the ouputs of Yk equations
    b=b+Yk;%using the dummy variable to be used in further summations 
end

Y=Y0+b;%now summing together the outputs of Y0 and Yk to be plotted on a graph
plot(x,Y,'r');%plotting the deflection (y axis) against distance along the beam (x-axis)
%giving relevent titles to the graph and the x and y axis
title(['Deflection due to ' int2str(num_of_loads) ' loads']);
xlabel('delfection (m)');
ylabel('distance from left hand end (m)');
grid on;