%Clear the command window.
clc;
%Clear all previous variables.
clear all;
%Close all previously opened figures or images.
close all;

%Initialize n, the size of the nonlinear system.
n = 2;

%Initialize some variables.
xcurrent = zeros(n,1);
F = zeros(n,1);
dF = zeros(n,1);
J = zeros(n,n);
deltax = zeros(n,1);

%Open an output file.
OutputFile = fopen('bailey_pset6_problem1_a_test.txt','w');

%Print title to the screen and to the output file.
fprintf('\n OUTPUT FROM bailey_pset6_problem1_a_test.m \n\n');
fprintf(OutputFile, '\n OUTPUT FROM bailey_pset6_problem1_a_test.m \n\n');

%Initialize some more variables.
it = 0;                   %Iteration counter (outer loop)
kmax = 25;               %Maximum number of iterations (outer loop)
epsilon = 0.00001;       %Tolerance
check = 10*epsilon;    %To be used in the outer While loop
c = 10^(-8);


%Ask user for initial guess vector.
xcurrent(1) = input(' Please input an initial guess for x1. \n');
xcurrent(2) = input(' Please input an initial guess for x2. \n');

%Print information about the method and the problem to the screen and to the output file.
fprintf('  Solving 2 Nonlinear Equations Using Generalized Newton Method (Numerical) \n');
fprintf('     Initial Guesses: x1 = %+1.5e and x2 = %+1.5e \n\n', xcurrent(1), xcurrent(2));
fprintf(OutputFile, '  Solving 2 Nonlinear Equations Using Generalized Newton Method (Numerical) \n');
fprintf(OutputFile, '     Initial Guesses: x1 = %+1.5e and x2 = %+1.5e \n\n', xcurrent(1), xcurrent(2));

%Print the column headings for the results table.
fprintf('%8s%4s%6s%10s%15s%16s\n', 'Iter.', 'x1', 'x2', 'check1');
fprintf(OutputFile, '%8s%4s%6s%10s%15s%16s\n', 'Iter.', 'x1', 'x2', 'check1');

%Print a horizontal line below the column headings.
fprintf('%s\n','------------------------------------------------------------------');
fprintf(OutputFile, '%s\n','------------------------------------------------------------------');

%Outer iteration loop (Nonlinear Jacobi)
while (check >= epsilon) && (it < kmax)
    
    %residual
    F(1) = f1(xcurrent(1),xcurrent(2));
    F(2) = f2(xcurrent(1),xcurrent(2));
    
    %jacobian
    theta = c*abs(xcurrent(1)) + c;
    dF(1) = f1(xcurrent(1)+ theta ,xcurrent(2));
    dF(2) = f2(xcurrent(1)+ theta ,xcurrent(2));
    for i = 1:n
        J(i,1) = (dF(i) - F(i))/theta;
    end
    
    %Recalculate theta for theta2
    theta = c*abs(xcurrent(2)) + c;
    dF(1) = f1(xcurrent(1),xcurrent(2) + theta);
    dF(2) = f2(xcurrent(1),xcurrent(2) + theta);
    for i = 1:n
        J(i,2) = (dF(i) - F(i))/theta;
    end
    
    %Using Gaussian Elimination to solve J*deltax =(-F);
    m = J(2,1)/J(1,1);
    J(2,2) = J(2,2) - (m * J(1,2));
    F(2) = F(2) - (m * F(1));
    
    %Calculate the Back Solve
    deltax(2) = -F(2)/J(2,2);
    deltax(1) = (-F(1)-(J(1,2)*deltax(2)))/J(1,1);
    
    xold = xcurrent;
    for i=1:n
        xcurrent(i) = xcurrent(i) + deltax(i);
    end
    it = it + 1;
    
    %Update check here
    maxtop = abs(xcurrent(1) - xold(1));
    for i = 1:n-1
        if abs(xcurrent(i+1) - xold(i+1)) > maxtop
            maxtop = abs(xcurrent(i+1) - xold(i+1));
        end
    end
    
    maxbot = abs(xcurrent(1));
    for i = 1:n-1
        if abs(xcurrent(i+1)) > maxbot
            maxbot = abs(xcurrent(i+1));
        end
    end
    check = maxtop/maxbot;
    
    %Print info about this iteration of Nonlinear Jacobi.
    fprintf('    %2d   %+1.5e   %+1.5e   %+1.5e\n', it, xcurrent(1), xcurrent(2), check);
    fprintf(OutputFile, '    %2d   %+1.5e   %+1.5e   %+1.5e\n', it, xcurrent(1), xcurrent(2), check);
    %Update variables in preparation for next iteration.
    %You may want to update additional variables here...
    
end

if check < epsilon
    fprintf('x is the converged solution to F(x) = 0 in %2d Newton iterations.',it);
    fprintf(OutputFile, 'x is the converged solution to F(x) = 0 in %2d Newton iterations.',it);
else
    fprintf('Converged solution not found within kmax iterations.');
    fprintf(OutputFile, 'Converged solution not found within kmax iterations.');
end

%Close the output file.
fclose(OutputFile);
 
%Define the functions
function answer = f1(x1,x2)

answer = (sin(x1*x2)/2)-(x1/2)-(x2/(4*pi));

return
end
function answer = f2(x1,x2)

answer = ((1-(1/(4*pi))) * (exp(2*x1) - exp(1))) - (((2*x1)-(x2/pi))*exp(1));

return
end
