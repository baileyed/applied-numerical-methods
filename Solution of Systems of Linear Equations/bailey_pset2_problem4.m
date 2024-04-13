function bailey_pset2_problem4

%Clear the command window.
clc;
%Clear all previous variables.
clear all;
%Close all previously opened figures or images.
close all;

%Initialize n, the size of the linear system.
n = 8;

%Open an output file.
OutputFile = fopen('bailey_pset2_problem4.txt','w');

%Print title to the screen and to the output file.
fprintf('\n OUTPUT FROM bailey_pset2_problem4.m \n\n');
fprintf(OutputFile, '\n OUTPUT FROM bailey_pset2_problem4.m \n\n');

%Choose which system to solve.
prompt = ' Please enter ''1'' to solve Afrac * x = bfrac. \n Please enter ''2'' to solve Aint * x = bint. \n';
RunChoice = input(prompt); 

%Print information about the method and the problem to the screen and to the output file.
%Initialize coefficient matrix A and RHS vector b.
if RunChoice == 1
    fprintf(' Solving Afrac * x = bfrac using Gaussian Elimination \n');
    fprintf(OutputFile, ' Solving Afrac * x = bfrac using Gaussian Elimination \n');
    [A, b] = GetFracSystem(n);
else
    fprintf(' Solving Aint * x = bint using Gaussian Elimination \n');
    fprintf(OutputFile, ' Solving Aint * x = bint using Gaussian Elimination \n');
    [A, b] = GetIntSystem(n);
end

%Print the current version of A and b to the screen and to the output file.
PrintAb(OutputFile, n, A, b);

for k = 1:(n-1)
    for i = (k+1):n
        if A(i,k)> A(k,k)
            Aik = A(i,k);
            A(i,k) = A(k,k);
            A(k,k) = Aik;
            fprintf(' \n The values %10.5f and %10.5f in positions (%2d, %2d) and (%2d, %2d) were switched. \n', A(i,k), A(k,k),i,k,k,k);
            fprintf(OutputFile, ' \n The values %10.5f and %10.5f in positions (%2d, %2d) and (%2d, %2d) were switched. \n', A(i,k), A(k,k),i,k,k,k);
        else
            fprintf(' \n No rows were swapped. \n');
            fprintf(OutputFile, ' \n No rows were swapped. \n');
        end
        m_ik = A(i,k)/A(k,k);
        A(i,k) = 0;
        for j = (k+1):n
            A(i,j) = A(i,j) - (m_ik * A(k,j));
        end
        b(i) = b(i) - (m_ik * b(k));
    end
    %print the updated A and the updated b
    PrintAb(OutputFile, n, A, b);
end

%Calculate the Back Solve
for i = n:-1:1
   x(i)= b(i)/A(i,i);
   for j = (k+1):n
       x(k) = (b(k) - A(k,j)*x(j))/A(k,k);
   end
end

for i = 1:n
   %Print the solution vector x.
   fprintf('\n x%2d = %10.5f \n', i, x(i));
   fprintf(OutputFile, '\n x%2d = %10.5f \n', i, x(i));
end

%Close the output file.
fclose(OutputFile);
end

function [Afrac, bfrac] = GetFracSystem(n)

Afrac = zeros(n,n);
bfrac = zeros(n,1);

Afrac(1,1) = -19/12;
Afrac(1,2) = 0.5;
Afrac(1,8) = 1/12;

Afrac(2,1) = 0.5;
Afrac(2,2) = -61/66;
Afrac(2,3) = 1/3;
Afrac(2,7) = 1/11;

Afrac(3,2) = 1/3;
Afrac(3,3) = -41/60;
Afrac(3,4) = 0.25;
Afrac(3,6) = 0.1;

Afrac(4,3) = 0.25;
Afrac(4,4) = -0.45;
Afrac(4,5) = 0.2;

Afrac(5,4) = 0.2;
Afrac(5,5) = -11/30;
Afrac(5,6) = 1/6;

Afrac(6,3) = 0.1;
Afrac(6,5) = 1/6;
Afrac(6,6) = -43/105;
Afrac(6,7) = 1/7;

Afrac(7,2) = 1/11;
Afrac(7,6) = 1/7;
Afrac(7,7) = -221/616;
Afrac(7,8) = 0.125;

Afrac(8,1) = 1/12;
Afrac(8,7) = 0.125;
Afrac(8,8) = -23/72;

bfrac(1) = -200;

return
end

function [Aint, bint] = GetIntSystem(n)

Aint = zeros(n,n);
bint = zeros(n,1);

Aint(1,1) = -19;
Aint(1,2) = 6;
Aint(1,8) = 1;

Aint(2,1) = 33;
Aint(2,2) = -61;
Aint(2,3) = 22;
Aint(2,7) = 6;

Aint(3,2) = 20;
Aint(3,3) = -41;
Aint(3,4) = 15;
Aint(3,6) = 6;

Aint(4,3) = 5;
Aint(4,4) = -9;
Aint(4,5) = 4;

Aint(5,4) = 6;
Aint(5,5) = -11;
Aint(5,6) = 5;

Aint(6,3) = 21;
Aint(6,5) = 35;
Aint(6,6) = -86;
Aint(6,7) = 30;

Aint(7,2) = 56;
Aint(7,6) = 88;
Aint(7,7) = -221;
Aint(7,8) = 77;

Aint(8,1) = 6;
Aint(8,7) = 9;
Aint(8,8) = -23;

bint(1) = -2400;
return
end


function PrintAb(OutputFile, n, A, b)

fprintf('\n Current version of A and b:\n\n');
fprintf(OutputFile, '\n Current version of A and b:\n\n');

for i = 1:n
    for j = 1:n
        fprintf(' %10.5f', A(i,j));
        fprintf(OutputFile, ' %10.5f', A(i,j));
    end
    fprintf(' | %10.5f\n', b(i));
    fprintf(OutputFile, ' | %10.5f\n', b(i));
end
fprintf(' \n');
fprintf(OutputFile, ' \n');

return
end
