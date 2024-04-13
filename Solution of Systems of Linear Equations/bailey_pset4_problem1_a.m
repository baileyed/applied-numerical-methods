function bailey_pset4_problem1_a

%Clear the command window.
clc;
%Clear all previous variables.
clear all;
%Close all previously opened figures or images.
close all;

%Initialize n, the size of the linear system.
n = 20;

%Open an output file.
OutputFile = fopen('bailey_pset4_problem1_a.txt','w');

%Print title to the screen and to the output file.
fprintf('\n OUTPUT FROM bailey_pset4_problem1_a.m \n\n');
fprintf(OutputFile, '\n OUTPUT FROM bailey_pset4_problem1_a.m \n\n');

%Print information about the method and the problem to the screen and to the output file.
fprintf(' Solving Z * x = b for a tridiagonal system of size n = %4d using Jacobi Iteration \n', n);
fprintf(OutputFile, ' Solving Z * x = b for a tridiagonal system of size n = %4d using Jacobi Iteration \n', n);

%Initialize the tridiagonal coefficient matrix Z and the RHS vector b.
[Zsub, Zmain, Zsuper, b] = GetZb(n);

%Print the column headings for the results table.
fprintf('\n k+1     x(1)         x(2)         x(3)        x(4)          x(n-3)       x(n-2)       x(n-1)       x(n)       check       alpha\n');
fprintf(OutputFile, '\n k+1     x(1)         x(2)         x(3)        x(4)          x(n-3)       x(n-2)       x(n-1)       x(n)       check       alpha\n');

%Print a horizontal line below the column headings.
fprintf('%s\n','----------------------------------------------------------------------------------------------------------------------------------');
fprintf(OutputFile, '%s\n','----------------------------------------------------------------------------------------------------------------------------------');

xnew = zeros(n,1);
xold = b;
epsilon = .0001;
check = 1;
k = 0;
kmax = 2500;
x = ones(n,1)*(1/(n+1)^2);

%Main Loop
while k < kmax && check >= epsilon
    if k+1 > kmax
        fprintf(' \n The maximum k value was reached. \n');
        fprintf(OutputFile, ' \n The maximum k value was reached. \n');
        break
    end
    for i=1:n
        if i == 1
            xnew(i) = ((b(i) - Zsuper(i)*xold(i+1))/Zmain(i));
        elseif i == n
            xnew(i) = ((b(i) - Zsub(i)*xold(i-1))/Zmain(i));
        else
            xnew(i) = ((b(i) - Zsub(i)*xold(i-1) - Zsuper(i)*xold(i+1))/Zmain(i));
        end
    end
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
    
    check = maxtop/maxbot;    
    sumnew = 0;
    sumold = 0;
    
    for i = 1:n
        sumnew = sumnew + ((xnew(i) - x(i))^2);
        sumold = sumold + ((xold(i) - x(i))^2);
        alpha = sqrt(sumnew)/sqrt(sumold);    
    end
    
    k = k + 1;    
    fprintf(' %2d ', k);
    fprintf(OutputFile, ' %2d ', k);
    
    for i=1:4
        fprintf(' %+1.4e ', xnew(i));
        fprintf(OutputFile, ' %+1.4e ', xnew(i));
    end
    
    for i=n-3:n
        fprintf(' %+1.4e ', xnew(i));
        fprintf(OutputFile, ' %+1.4e ', xnew(i));    
    end
    
    fprintf('%+1.4e %+1.4e\n', check, alpha);
    fprintf(OutputFile, '%+1.4e %+1.4e\n', check, alpha);
    
    for i=1:n
        xold(i) = xnew(i);
    end
end

%Print a conclusion statement.
if check >= epsilon
    fprintf(' The method did not converge in %2d iterations.\n\n', k);
    fprintf(OutputFile, ' The method did not converge in %2d iterations.\n\n', k);
else
    fprintf(' A converged solution was found after %2d iterations.\n\n', k);
    fprintf(OutputFile, ' A converged solution was found after %2d iterations.\n\n', k);
end

%Close the output file.
fclose(OutputFile);

end
function [Zsub, Zmain, Zsuper, b] = GetZb(n)

%Initialize commonly used variables to avoid unnecessarily repetitive
%computations when initializing the coefficient matrix.
np1sq = (n+1)*(n+1);
tempmain = 2*np1sq;  %Value on main diagonal of Z
tempoff = -np1sq;    %Value on off-diagonals of Z

%Initialize the tridiagonal coefficient matrix Z.  Notice that it is possible
%to assign all entries in Z without using any conditionals (any "if"s).
%Avoiding unnecessary conditionals (especially inside of loops, where
%they might execute MANY times) is a good programming habit.
Zmain = zeros(n,1);
Zsuper = zeros(n,1);
Zsub = zeros(n,1);

Zmain(1) = tempmain;
Zsuper(1) = tempoff;
for i=2:n-1
    Zsub(i) = tempoff;
    Zmain(i) = tempmain;
    Zsuper(i) = tempoff;
end
Zsub(n) = tempoff;
Zmain(n) = tempmain;

%Initialize the RHS vector b.
b = zeros(n,1);
b(1) = 1;
b(n) = 1;

return
end
