function bailey_pset1_c

%Clear the command window.
clc;
%Clear all previous variables.
clear all;
%Close all previously opened figures or images.
close all;

%Open an output file.
OutputFile = fopen('bailey_pset1_c.txt','w');

%Print title to the screen and to the output file.
fprintf('\n OUTPUT FROM bailey_pset1_c.m \n\n');
fprintf(OutputFile, '\n OUTPUT FROM bailey_pset1_c.m \n\n');

function answer = f(x)

answer = x^2-4*x*(1 + cos(2*x)) + 2;

return
end

%Initialize endpoints of interval.
x_1 = input(' Please input the left endpoint, x_1. \n');   %Left endpoint of interval
x_0 = input(' Please input the right endpoint, x_0. \n');  %Right endpoint of interval

%Initialize some other variables.
n = 0;                %Iteration counter
nmax = 50;            %Maximum allowed number of iterations
epsilon = 0.00001;    %Tolerance
check = 10 * epsilon; %To be used in the While loop

fx_1 = f(x_1);
fx_0 = f(x_0);
if fx_1*fx_0 > 0
    fprintf('There is a root \n');
    exit()
end

%Print information about the method and the problem to the screen and to the output file.
fprintf(' Solving a nonlinear equation using Secant Method \n');
fprintf(' with initial guesses x_1 = %9.6f and x_0 = %9.6f \n\n', x_1, x_0);
fprintf(OutputFile, ' Solving a nonlinear equation using Secant Method \n');
fprintf(OutputFile, ' with initial guesses x_1 = %9.6f and x_0 = %9.6f \n\n', x_1, x_0);

%Print the column headings for the results table.
fprintf('%10s%22s%28s\n', 'Iteration', 'x', 'f(x)');
fprintf(OutputFile, '%10s%22s%28s\n', 'Iteration', 'x', 'f(x)');

%Print a horizontal line below the column headings.
fprintf('%s\n','------------------------------------------------------------------');
fprintf(OutputFile, '%s\n','------------------------------------------------------------------');

%Main loop
while (check >= epsilon) && (n < nmax)
    if n+1 > nmax
        fprintf('The nmax has been reached \n')
        exit()
    end
    %We set xnew to a random number in the interval [0,1].
    xnew = x_0-(fx_0/((fx_0-fx_1)/(x_0-x_1)));
    fxnew = f(xnew);
    check = abs(fxnew);
    %Now we print the iteration number n+1, the value of xnew, and
    %the value of check (the quantity that shows whether we're making
    %progress toward meeting the termination criterion).
    fprintf('    %2d                  %12.6f               %12.6f\n', n+1, xnew, fxnew);
    fprintf(OutputFile, '    %2d                  %12.6f               %12.6f\n', n+1, xnew, fxnew);
    x_1 = x_0;
    fx_1 = fx_0;
    x_0 = xnew;
    fx_0 = fxnew;
    n = n+1;
end

%Print another horizontal line.
fprintf('%s\n','------------------------------------------------------------------');
fprintf(OutputFile, '%s\n','------------------------------------------------------------------');

%Print a conclusion statement.
if check >= epsilon
    fprintf(' The method did not converge in %2d iterations.\n\n', n);
    fprintf(OutputFile, ' The method did not converge in %2d iterations.\n\n', n);
else
    fprintf(' The method converged to %12.6f after %2d iterations.\n\n', xnew, n);
    fprintf(OutputFile, ' The method converged to %12.6f after %2d iterations.\n\n', xnew, n);
end

%Close the output file.
fclose(OutputFile);
end