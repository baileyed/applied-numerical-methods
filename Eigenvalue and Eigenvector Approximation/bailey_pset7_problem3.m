function bailey_pset7_problem3

%Clear the command window.
clc;
%Clear all previous variables.
clear all;
%Close all previously opened figures or images.
close all;

%Initialize n, the size of the linear system.
n = 20;

%Open an output file.
OutputFile = fopen('bailey_pset7_problem3.txt','w');

%Print title to the screen and to the output file.
fprintf('\n OUTPUT FROM bailey_pset7_problem3.m \n\n');
fprintf(OutputFile, '\n OUTPUT FROM bailey_pset7_problem3.m \n\n');

%Print information about the method and the problem to the screen and to the output file.
fprintf(' Finding all the eigenvalues and eigenvectors of A for a system of size n = %4d using Jacobi’s Method \n', n);
fprintf(OutputFile, ' Finding all the eigenvalues and eigenvectors of A for a system of size n = %4d using Jacobi’s Method \n', n);

%Initialize the dense coefficient matrix A and the RHS vector b.
[A] = GetAb(n);
k = 0;
kmax = 100;
epsilon = .00001;

%Print the column headings for the results table.
fprintf('\n k +1      x(1)         x(2)         x(3)         x(4)        x(5)       x(6)       x(7)        x(8)           x(9)        x(10)       x(11)         x(12)         x(13)         x(14)        x(15)       x(16)       x(17)        x(18)         x(19)          x(20)       offA\n');
fprintf(OutputFile, '\n k +1      x(1)         x(2)         x(3)         x(4)        x(5)       x(6)       x(7)        x(8)           x(9)        x(10)       x(11)         x(12)         x(13)         x(14)        x(15)       x(16)       x(17)        x(18)         x(19)          x(20)       offA\n');

%Print a horizontal line below the column headings.
fprintf('%s\n','-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');
fprintf(OutputFile, '%s\n','-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------');

E = zeros(n,n);
for i = 1:n
        E(i,i) = 1;
end

offA = 0;
for i = 1:n
    for j = 1:n
        offA = offA + (A(i,j)^2);
    end
    offA = offA - (A(i,i)^2);
end


while (offA >= epsilon) && (k < kmax)
    
    p = 0;
    q = 0;
    maxA = 0;
    for i = 1:n
        for j = 1:n
            if (maxA < abs(A(i,j))) && (i<j)
                maxA = abs(A(i,j));
                p = i;
                q = j;
            end
        end
    end
    
    apq = A(p,q);
    
    tau = (A(q,q) - A(p,p))/(2*A(p,q));
    if tau >= 0
        signt = 1;
    else
        signt = -1;
    end
    t = signt/(abs(tau) + sqrt(1 + (tau)^2));
    c = 1/(sqrt(1 + t^2));
    s = t*c;
    
    A(p,p) = A(p,p) - t*A(p,q);
    A(q,q) = A(q,q) + t*A(p,q);
    
    for j = 1:(p-1)
        ajp = A(j,p);
        A(j,p) = c*ajp - s*A(j,q);
        A(j,q) = s*ajp + c*A(j,q);
    end
    for j = (p+1):(q-1)
        apj = A(p,j);
        A(p,j) = c*apj - s*A(j,q);
        A(j,q) = s*apj + c*A(j,q);
    end
    for j = (q+1):n
        apj = A(p,j);
        A(p,j) = c*apj - s*A(q,j);
        A(q,j) = s*apj + c*A(q,j);
    end
    
    A(p,q) = 0;
    
    for i = 1:n
        for j = 1:n
            if i>j
                A(i,j) = A(j,i);
            end
        end
    end
    
    for j = 1:n
        eip = E(j,p);
        E(j,p) = c*eip - s*E(j,q);
        E(j,q) = s*eip + c*E(j,q);
    end
    
    offAold = offA;
    offA = offA - 2*(apq)^2;
    
    alpha = offA/offAold;
    
    k = k + 1;
    
    %Print Statements
    
    fprintf('\n %2d',k);
    fprintf(OutputFile, '\n %2d',k);
    
    for i = 1:10
        fprintf('  %1.5e', A(i,i));
        fprintf(OutputFile, '  %1.5e', A(i,i));
    end
    for i = 11:20
        fprintf('  %1.5e', A(i,i));
        fprintf(OutputFile, '  %1.5e', A(i,i));
    end
    fprintf('  %1.5e', offA);
    fprintf(OutputFile, '  %1.5e', offA);
    
end

if offA < epsilon
    fprintf('\n Eigenvalue found in %2d iterations with corresponding eigenvector x(k)', k);
    fprintf(OutputFile, '\n Eigenvalue found in %2d iterations with corresponding eigenvector x(k)', k);
    fprintf('\n Eigenvalues of A appear along the diagonal of A(k). \n');
    fprintf(OutputFile, '\n Eigenvalues of A appear along the diagonal of A(k).\n');
    for i=1:10
        fprintf(' %+1.5e', A(i,i));
        fprintf(OutputFile, ' %+1.5e', A(i,i));
    end
    fprintf('\n');
    fprintf(OutputFile, '\n');
    for i=11:n
        fprintf(' %+1.5e', A(i,i));
        fprintf(OutputFile, ' %+1.5e', A(i,i));
    end
    fprintf('\n Eigenvectors of A are the columns of E(k). \n');
    fprintf(OutputFile, '\n Eigenvectors of A are the columns of E(k). \n');
    for i = 1:n
        fprintf(' v(%2d)',i);
        fprintf(OutputFile, ' v(%2d)',i);
        for j = 1:n
            evector(j) = E(j,i);
            fprintf(' %+1.5e', evector(j));
            fprintf(OutputFile, ' %+1.5e', evector(j));
        end
        fprintf(' \n');
        fprintf(OutputFile, ' \n');
    end
else
    fprintf('\n Eigenvalues and vectors not found within kmax iterations. End.');
    fprintf(OutputFile, '\n Eigenvalues and vectors not found within kmax iterations. End.');
end

%Close the output file.
fclose(OutputFile);

end

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