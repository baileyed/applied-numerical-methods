function bailey_pset6_problem3

%Clear the command window.
clc;
%Clear all previous variables.
clear all;
%Close all previously opened figures or images.
close all;

%Initialize n, the size of the linear system.
n = 20;

%Open an output file.
OutputFile = fopen('bailey_pset6_problem3.txt','w');

%Print title to the screen and to the output file.
fprintf('\n OUTPUT FROM bailey_pset6_problem3.m \n\n');
fprintf(OutputFile, '\n OUTPUT FROM bailey_pset6_problem3.m \n\n');

%Print information about the method and the problem to the screen and to the output file.
fprintf(' Solving A * x = b for a dense system of size n = %4d using Conjugate Gradient Method \n', n);
fprintf(OutputFile, ' Solving A * x = b for a dense system of size n = %4d using Conjugate Gradient Method \n', n);

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
kmax = n;
epsilon = .0001;

%initialize rold
rold = zeros(1,n);
for i = 1:n
    Ap = 0;
    for j = 1:n
        Ap = Ap + A(i,j)*xnew(j);
    end
    rold(i) = b(i) - Ap;
end

rcurr = rold;
rnorm = 0;
for i = 1:n
    rnorm =  rnorm + rcurr(i)^2;
end
rtransr = rnorm;
rnorm = sqrt(rnorm);

while k < kmax && rnorm >= epsilon
    %rtranspose*r for rold
    roldtr = 0;
    for i = 1:n
        roldtr =  roldtr + rold(i)^2;
    end
    
    %Find delta and pnew
    if k == 0
        delta = 0;
        pnew = rcurr;
    else
        delta = (rtransr)/(roldtr);
        for i = 1:n
            pnew(i) = rcurr(i) + delta*pnew(i);
        end
    end
    
    %Calculate
    ptransA = 0;
    for i = 1:n
        for j = 1:n
            ptransA = ptransA + A(i,j)*pnew(i)*pnew(j);
        end
    end
    
    gamma = (rtransr)/(ptransA);
    
    for i = 1:n
        xnew(i) = xnew(i) + gamma*pnew(i);
    end
    
    rold = rcurr;
    
    for i = 1:n
        Ap = 0;
        for j = 1:n
            Ap = Ap + A(i,j)*pnew(j);
        end
        rcurr(i) = rold(i) - gamma*Ap;
    end
    
    k = k+1;
    
    rnorm = 0;
    for i = 1:n
        rnorm =  rnorm + rcurr(i)^2;
    end
    rtransr = rnorm;
    rnorm = sqrt(rnorm);
    
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
    
    fprintf(' %+1.4e \n', rnorm);
    fprintf(OutputFile, ' %+1.4e \n', rnorm);
end

if rnorm < epsilon
    fprintf('x is the converged solution to Ax = b in %2d iterations. \n',k);
    fprintf(OutputFile, 'x is the converged solution to Ax = b in %2d iterations. \n',k);
else
    fprintf('Converged solution not found within kmax iterations. \n');
    fprintf(OutputFile, 'Converged solution not found within kmax iterations. \n');
end

%Print the column headings for the circle table.
fprintf('\n s      center         radius          ratio\n');
fprintf(OutputFile, '\n s      center         radius          ratio\n');

%Print a horizontal line below the column headings.
fprintf('%s\n','------------------------------------------------------------------------------');
fprintf(OutputFile, '%s\n','------------------------------------------------------------------------------');

for i = 1:n
    radius = 0;
    for j = 1:n
        if i ~= j
        radius = radius + abs(A(i,j));
        end
    end
    center = A(i,i);
    ratio = radius/center;
    fprintf('\n %2d %+1.8e %+1.8e %+1.8e', i, center, radius, ratio);
    fprintf(OutputFile,'\n %2d %+1.8e %+1.8e %+1.8e', i, center, radius, ratio);
end

%Print the column headings for the eigen table.
fprintf('\n s       lambda          hsame           hdiff               v(1)          v(2)         v(3)         v(4)        v(n-3)        v(n-2)       v(n-1)       v(n)\n');
fprintf(OutputFile, '\n s       lambda          hsame           hdiff               v(1)          v(2)         v(3)         v(4)        v(n-3)        v(n-2)       v(n-1)       v(n)\n');

%Print a horizontal line below the column headings.
fprintf('%s\n','--------------------------------------------------------------------------------------------------------------------------------------------------------------');
fprintf(OutputFile, '%s\n','--------------------------------------------------------------------------------------------------------------------------------------------------------------');

vector = zeros(1,n);

xsame = zeros(1,n);
xsame(1) = 1;
xsame(n) = 1;

xdiff = zeros(1,n);
xdiff(1) = 1;
xdiff(n) = -1;
for s = 1:n
    lambda = 1/(((n+1)^2)*(2 - 2*cos((s*pi)/(n+1))));
    
    for t = 1:n
        vector(t) = sin((s*t*pi)/(n+1));
    end
    
    d = 0;
    diffn = 0;
    samen = 0;
    for i = 1:n
        for j = 1:n
            d = d + vector(i)*vector(j);
            samen = samen + vector(i)*xsame(j);
            diffn = diffn + vector(i)*xdiff(j);
        end
        
    end
    
    hsame = samen/d;
    hdiff = diffn/d;
%     samen = 0;
%     for i = 1:n
%         for j = 1:n
%             samen = samen + vector(i)*xsame(j);
%         end
%     end
%     
%     diffn = 0;
%     for i = 1:n
%         for j = 1:n
%             diffn = diffn + vector(i)*xdiff(j);
%         end
%     end
    
    
    
    fprintf('\n %2d %+1.8e', s, lambda);
    fprintf(OutputFile,'\n %2d  %+1.8e  ', s, lambda);
%     
%     for t=1:4
%         fprintf(' %12.5f', vector(t));
%         fprintf(OutputFile, ' %12.5f', vector(t));
%     end
%     for t=n-3:n
%         fprintf(' %12.5f', vector(t));
%         fprintf(OutputFile, ' %12.5f', vector(t));
%     end
% 
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