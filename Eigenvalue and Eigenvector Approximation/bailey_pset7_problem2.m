%Clear the command window.
clc;
%Clear all previous variables.
clear all;
%Close all previously opened figures or images.
close all;

%Initialize n, the size of the nonlinear system.
n = 20;

%Initialize some variables.
[A] = GetAb(n);

%Open an output file.
OutputFile = fopen('bailey_pset7_problem2.txt','w');

%Print title to the screen and to the output file.
fprintf('\n OUTPUT FROM bailey_pset7_problem2.m \n\n');
fprintf(OutputFile, '\n OUTPUT FROM bailey_pset7_problem2.m \n\n');

%Initialize some more variables.
k = 0;                   %Iteration counter (outer loop)
kmax = 8000;               %Maximum number of iterations (outer loop)
epsilon = 0.00001;       %Tolerance

%Ask user for initial guess vector.
prompt = ' Please enter ''1'' to use x same. \n Please enter ''2'' to use x different. \n';
RunChoice = input(prompt);

%Print information about the method and the problem to the screen and to the output file.
fprintf(' Finding Eigenvectors and Eigenvalues Using Shifted Power Method \n');
fprintf(OutputFile, ' Finding Eigenvectors and Eigenvalues Using Shifted Power Method \n');

%Print information about the method and the problem to the screen and to the output file.
%Initialize coefficient matrix A and RHS vector b.
if RunChoice == 1
    fprintf(' Using x same as the initial guess vector. \n');
    fprintf(OutputFile, ' Using x same as the initial guess vector. \n');
    xsame = zeros(1,n);
    xsame(1) = 1;
    xsame(n) = 1;
    xcurrent = xsame;
else
    fprintf(' Using x different as the initial guess vector. \n');
    fprintf(OutputFile, ' Using x different as the initial guess vector. \n');
    xdiff = zeros(1,n);
    xdiff(1) = 1;
    xdiff(n) = -1;
    xcurrent = xdiff;
end

%Print the column headings for the results table.
fprintf('%4s%8s%10s%14s%14s%14s%12s%12s%12s%14s%12s\n', 'Iter.', 'x(1)', 'x(2)', 'x(3)', 'x(4)', 'x(n-3)', 'x(n-2)', 'x(n-1)', 'x(n)', 'lambda', 'check');
fprintf(OutputFile, '%4s%8s%10s%14s%14s%14s%12s%12s%12s%14s%12s\n', 'Iter.', 'x(1)', 'x(2)', 'x(3)', 'x(4)', 'x(n-3)', 'x(n-2)', 'x(n-1)', 'x(n)', 'lambda', 'check');

%Print a horizontal line below the column headings.
fprintf('%s\n','------------------------------------------------------------------------------------------------------------------------------------');
fprintf(OutputFile, '%s\n','------------------------------------------------------------------------------------------------------------------------------------');

%compute inital values
sigma = 1/(((n+1)^2)*(2-2*cos(pi/(n+1))));
M = A;
for i = 1:n
    M(i,i) = A(i,i) - sigma;
end

%divide by 2-norm to get unit vector
xnorm = 0;
for i = 1:n
    xnorm =  xnorm + abs(xcurrent(i))^2;
end
xunit = sqrt(xnorm);

for i = 1:n
    xcurrent(i) = xcurrent(i)/xunit;
end

%pcurrent = xcurrent*A --> v-m mult
pcurrent = zeros(n,1);
for i = 1:n
    Mx = 0;
    for j = 1:n
        Mx = Mx + M(i,j)*xcurrent(j);
    end
    pcurrent(i) = Mx;
end

%lambda = xcurrent*pcurrent --> dot product
lambda = 0;
for i = 1:n
    lambda = lambda + xcurrent(i)*pcurrent(i);
end

%residual = lambda*xcurrent - pcurrent
residual = 0;
rnorm = 0;
for i = 1:n
    residual = lambda*xcurrent(i) - pcurrent(i);
    rnorm =  rnorm + abs(residual)^2;
end
check = sqrt(rnorm);

while (check >= epsilon) && (k < kmax)
    
    %divide by 2-norm to get unit vector
    pnorm = 0;
    for i = 1:n
        pnorm =  pnorm + abs(pcurrent(i))^2;
    end
    punit = sqrt(pnorm);
    
    for i = 1:n
        xcurrent(i) = pcurrent(i)/punit;
    end
    
    %pcurrent = xcurrent*A --> v-m mult
    for i = 1:n
        Mx = 0;
        for j = 1:n
            Mx = Mx + M(i,j)*xcurrent(j);
        end
        pcurrent(i) = Mx;
    end
    
    %lambda = xcurrent*pcurrent --> dot product
    lambda = 0;
    for i = 1:n
        lambda = lambda + xcurrent(i)*pcurrent(i);
    end
    
    %residual = lambda*xcurrent - pcurrent    
    %Update check here - the 2 norm of the residual
    rnorm = 0;
    for i = 1:n
        residual = lambda*xcurrent(i) - pcurrent(i);
        rnorm =  rnorm + abs(residual)^2;
    end
    check = sqrt(rnorm);
    
    k = k + 1;
    
    %Print info about this iteration of Nonlinear Jacobi.
    fprintf('%2d', k);
    fprintf(OutputFile, '%2d', k);
    
    for i=1:4
        fprintf(' %+1.5e', xcurrent(i));
        fprintf(OutputFile, ' %+1.5e', xcurrent(i));
    end
    for i=n-3:n
        fprintf(' %+1.5e', xcurrent(i));
        fprintf(OutputFile, ' %+1.5e', xcurrent(i));
    end
    
    fprintf(' %+1.5e %+1.5e \n', (lambda + sigma), check);
    fprintf(OutputFile, ' %+1.5e %+1.5e \n', (lambda + sigma), check);
    
end

if check < epsilon
    fprintf(' Eigenvalue we seek, found in %2d iterations, is lambda + sigma = %+1.5e with corresponding eigenvector x(k)', k, (lambda + sigma));
    fprintf(OutputFile, ' Eigenvalue we seek, found in %2d iterations, is lambda + sigma = %+1.5e with corresponding eigenvector x(k)', k, (lambda + sigma));
    fprintf('\n All of the entries of the final eigenvector:\n');
    fprintf(OutputFile, '\n All of the entries of the final eigenvector:\n');
    for i=1:10
        fprintf(' %+1.5e', xcurrent(i));
        fprintf(OutputFile, ' %+1.5e', xcurrent(i));
    end
    fprintf('\n');
    fprintf(OutputFile, '\n');
    for i=11:n
        fprintf(' %+1.5e', xcurrent(i));
        fprintf(OutputFile, ' %+1.5e', xcurrent(i));
    end
else
    fprintf(' Converged solution not found within kmax iterations. End.');
    fprintf(OutputFile, ' Converged solution not found within kmax iterations. End.');
end

%Close the output file.
fclose(OutputFile);

function [A] = GetAb(n)

%Initialize a commonly used variable to avoid unnecessarily repetitive
%computations when initializing the coefficient matrix.
factor = 1./((n+1)^3);

%Initialize the dense coefficient matrix A.
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
return
end