%Clear the command window.
clc;
%Clear all previous variables.
clear all;
%Close all previously opened figures or images.
close all;

%Initialize n, the size of the nonlinear system.
n = 2;

%Initialize some variables.
xold = zeros(n,1);
xnew = zeros(n,1);

%Open an output file.
OutputFile = fopen('bailey_pset5_problem2_a.txt','w');

%Print title to the screen and to the output file.
fprintf('\n OUTPUT FROM bailey_pset5_problem2_a.m \n\n');
fprintf(OutputFile, '\n OUTPUT FROM bailey_pset5_problem2_a.m \n\n');

%Initialize some more variables.
k = 0;                   %Iteration counter (outer loop)
kmax = 25;               %Maximum number of iterations (outer loop)
lmax = 100;               %You should set an 'lmax' value, too.
epsilon = 0.00001;       %Tolerance
epsilon2 = 0.000001;
check1 = 10*epsilon;    %To be used in the outer While loop

ltot = 0;
%Ask user for initial guess vector.
xold(1) = input(' Please input an initial guess for x1. \n');
xold(2) = input(' Please input an initial guess for x2. \n');

%Print information about the method and the problem to the screen and to the output file.
fprintf('  Solving 2 Nonlinear Equations Using Nonlinear Jacobi Iteration \n');
fprintf('     Initial Guesses: x1 = %+1.5e and x2 = %+1.5e \n\n', xold(1), xold(2));
fprintf(OutputFile, '  Solving 2 Nonlinear Equations Using Nonlinear Jacobi Iteration \n');
fprintf(OutputFile, '     Initial Guesses: x1 = %+1.5e and x2 = %+1.5e \n\n', xold(1), xold(2));

%Print the column headings for the results table.
fprintf('%8s%4s%6s%10s%15s%16s\n', 'Iter.', 'l1', 'l2', 'x1', 'x2', 'check1');
fprintf(OutputFile, '%8s%4s%6s%10s%15s%16s\n', 'Iter.', 'l1', 'l2', 'x1', 'x2', 'check1');

%Print a horizontal line below the column headings.
fprintf('%s\n','------------------------------------------------------------------');
fprintf(OutputFile, '%s\n','------------------------------------------------------------------');

%Outer iteration loop (Nonlinear Jacobi)
while (check1 >= epsilon) && (k < kmax)
    
    check2 = 10*epsilon2;
    fx1 = f1(xold(1),xold(2));
    fpx1 = fp1(xold(1),xold(2));
    xold_inner = xold(1);
    l1 = 0;
    while (check2 >= epsilon2) && (l1 < lmax)
        %Use single-eq. Newton's Method to solve 1st eq. for 1st unknown.
        %Iteration counter (inner loop)
        xnew_inner(1) = xold_inner - (fx1/fpx1);
        fx1 = f1(xnew_inner(1),xold(2));
        fpx1 = fp1(xnew_inner(1),xold(2));
        
        l1 = l1 + 1;
        check2 = (abs(xnew_inner(1) - xold_inner))/(abs(xnew_inner(1)));
        xold_inner = xnew_inner(1);
    end
    
    check2 = 10*epsilon2;
    fx2 = f2(xold(1),xold(2));
    fpx2 = fp2(xold(1),xold(2));
    xold_inner = xold(2);
    l2 = 0;
    %Use single-eq. Newton's Method to solve 2nd eq. for 2nd unknown.
    %Iteration counter (inner loop)
    while (check2 >= epsilon2) && (l2 < lmax)
        xnew_inner(2) = xold_inner - (fx2/fpx2);
        fx2 = f2(xold(1),xnew_inner(2));
        fpx2 = fp2(xold(1),xnew_inner(2));
        
        l2 = l2 + 1;
        check2 = (abs(xnew_inner(2) - xold_inner))/(abs(xnew_inner(2)));
        xold_inner = xnew_inner(2);
    end
    
    ltot = ltot + l1 + l2;
    k = k + 1;
    xnew = xnew_inner;
    
    %Update check here
    maxtop = abs(xnew(1) - xold(1));
    for i = 1:n-1
        if abs(xnew(i+1) - xold(i+1)) > maxtop
            maxtop = abs(xnew(i+1) - xold(i+1));
        end
    end
    
    maxbot = abs(xnew(1));
    for i = 1:n-1
        if abs(xnew(i+1)) > maxbot
            maxbot = abs(xnew(i+1));
        end
    end
    check1 = maxtop/maxbot;
    xold = xnew;
    
    %Print info about this iteration of Nonlinear Jacobi.
    fprintf('    %2d   %3d   %3d   %+1.5e   %+1.5e   %+1.5e\n', k, l1, l2, xnew(1), xnew(2), check1);
    fprintf(OutputFile, '    %2d   %3d   %3d   %+1.5e   %+1.5e   %+1.5e\n', k, l1, l2, xnew(1), xnew(2), check1);
    %Update variables in preparation for next iteration.
    %You may want to update additional variables here...
    
end

if check1 < epsilon
    fprintf('x is the converged solution to Ax = b in %2d Newton iterations.',ltot);
     fprintf(OutputFile, 'x is the converged solution to Ax = b in %2d Newton iterations.',ltot);
else
    fprintf('Converged solution not found within kmax iterations.');
    fprintf(OutputFile, 'Converged solution not found within kmax iterations.');
end

%Close the output file.
fclose(OutputFile);

%Define the functions and their derivatives.
function answer = f1(x1,x2)

answer = (sin(x1*x2)/2)-(x1/2)-(x2/(4*pi));

return
end
function answer = fp1(x1,x2)

answer = ((x2*cos(x1*x2))/2)-(1/2);

return
end
function answer = f2(x1,x2)

answer = ((1-(1/(4*pi))) * (exp(2*x1) - exp(1))) - (((2*x1)-(x2/pi))*exp(1));

return
end
function answer = fp2(x1,x2)

answer = (exp(1)/pi);

return
end
