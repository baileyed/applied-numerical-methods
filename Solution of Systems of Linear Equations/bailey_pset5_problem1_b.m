function bailey_pset5_problem1_b

%Clear the command window.
clc;
%Clear all previous variables.
clear all;
%Close all previously opened figures or images.
close all;

%Initialize n, the size of the linear system.
n = 200;

%Open an output file.
OutputFile = fopen('bailey_pset5_problem1_b.txt','w');

%Print title to the screen and to the output file.
fprintf('\n OUTPUT FROM bailey_pset5_problem1_b.m \n\n');
fprintf(OutputFile, '\n OUTPUT FROM bailey_pset5_problem1_b.m \n\n');

%Print information about the method and the problem to the screen and to the output file.
fprintf(' Solving A * x = b for a dense system of size n = %4d using Steepest Descent \n', n);
fprintf(OutputFile, ' Solving A * x = b for a dense system of size n = %4d using Steepest Descent \n', n);

%Initialize the dense coefficient matrix A and the RHS vector b.
[A, b] = GetAb(n);

%Print the column headings for the results table.
fprintf('\n k+1      x(1)         x(2)         x(3)         x(4)        x(n-3)       x(n-2)       x(n-1)        x(n)        rnorm\n');
fprintf(OutputFile, '\n k+1      x(1)         x(2)         x(3)         x(4)        x(n-3)       x(n-2)       x(n-1)        x(n)        rnorm\n');

%Print a horizontal line below the column headings.
fprintf('%s\n','------------------------------------------------------------------------------------------------------------------------');
fprintf(OutputFile, '%s\n','------------------------------------------------------------------------------------------------------------------------');

%Fill the xnew vector with random numbers in order demonstrate how to print it.
xnew = b;
k = 0;
kmax = 40000;
epsilon = .0001;

%initialize rnew
r = zeros(n,1);
for i = 1:n
    Ax = 0;
    for j = 1:n
        Ax = Ax + A(i,j)*xnew(j);
    end
    r(i) = b(i) - Ax;
end

rnorm = 10^5;

while k < kmax && rnorm >= epsilon
    rnorm = 0;
    for i = 1:n
        rnorm =  rnorm + r(i)^2;
    end
    rtransr = rnorm;
    rnorm = sqrt(rnorm);
    
    rtransA = 0;
    for i = 1:n
        for j = 1:n
            rtransA = rtransA + A(i,j)*r(i)*r(j);
        end
    end
    
    gamma = (rtransr)/(rtransA);
    
    for i = 1:n
        xnew(i) = xnew(i) + gamma*r(i);
    end
    
    for i = 1:n
        Ax = 0;
        for j = 1:n
            Ax = Ax + A(i,j)*xnew(j);
        end
        r(i) = b(i) - Ax;
    end
    
    k = k+1;
    fprintf('%2d', k);
    fprintf(OutputFile, '%2d', k);
    for i=1:4
        fprintf(' %12.5f', xnew(i));
        fprintf(OutputFile, ' %12.5f', xnew(i));
    end
    for i=n-3:n
        fprintf(' %12.5f', xnew(i));
        fprintf(OutputFile, ' %12.5f', xnew(i));
    end
    
    fprintf(' %+1.4e', rnorm);
    fprintf(OutputFile, ' %+1.4e', rnorm);
    
    fprintf('\n');
    fprintf(OutputFile, '\n');
end

if rnorm < epsilon
    fprintf('x is the converged solution to Ax = b in %2d iterations.',k);
     fprintf(OutputFile, 'x is the converged solution to Ax = b in %2d iterations.',k);
else
    fprintf('Converged solution not found within kmax iterations.');
    fprintf(OutputFile, 'Converged solution not found within kmax iterations.');
end

%Close the output file.
fclose(OutputFile);

end

function [A, b] = GetAb(n)

%Initialize a commonly used variable to avoid unnecessarily repetitive
%computations when initializing the coefficient matrix.
factor = 1./((n+1)^3);

%Initialize the dense coefficient matrix A.  Notice that it is possible
%to assign all entries in A without using any conditionals (any "if"s).
%Avoiding unnecessary conditionals (especially inside of loops, where
%they might execute MANY times) is a good programming habit.
A = zeros(n,n);
for i = 1:n
    %First compute the diagonal entry in Row i.
    A(i,i) = i*(n+1-i)*factor;
    for j = i+1:n
        %Compute entries in Row i that are to the right of the diagonal.
        A(i,j) = A(i,j-1) - i*factor;
        %Force A to be symmetric.
        A(j,i) = A(i,j);
    end
end

%Initialize the RHS vector b.
b = ones(n,1);

return
end