function bailey_pset10_problem1

% Each function is defined as a comment at the bottom. To run each, remove the comment.

%Clear the command window.
clc;
%Clear all previous variables.
clear all;
%Close all previously opened figures or images.
close all;

%Initialize kmax, the maximum number of iterations.
kmax = 10;

%Open an output file.
OutputFile = fopen('bailey_pset10_problem1.txt','w');

%Print title to the screen and to the output file.
fprintf('\n OUTPUT FROM bailey_pset10_problem1.m \n\n');
fprintf(OutputFile, '\n OUTPUT FROM bailey_pset10_problem1.m \n\n');

%Initialize integration limits.
a = input(' Please input the lower limit of integration, a. \n');
b = input(' Please input the upper limit of integration, b. \n');
Actual = input(' Please input the actual solution to the integral. \n');

%Initialize some other variables.
k = 0;                %Iteration counter
epsilon = 10^(-6);   %Tolerance
check = 10 * epsilon; %To be used in the While loop

%Print information about the method and the problem to the screen and to the output file.
fprintf(' Using Romberg Integration to approximate the value of a definite integral from a = %5.2f to b = %5.2f\n\n', a, b);
fprintf(OutputFile, ' Using Romberg Integration to approximate the value of a definite integral from a = %5.2f to b = %5.2f\n\n', a, b);

%Compute initial value of h.
h = b - a;

%Get the values of all of the Romberg coefficients.
C = GetCvalues(kmax);

%Print the column headings for the results table.
fprintf('%6s%8s%12s%12s%12s%12s%12s%11s\n', 'Iter.', 'Col. 0', 'Col. 1', 'Col. 2', 'Col. 3', 'Col. 4', 'Col. 5', '...');
fprintf(OutputFile, '%6s%8s%12s%12s%12s%12s%12s%11s\n', 'Iter.', 'Col. 0', 'Col. 1', 'Col. 2', 'Col. 3', 'Col. 4', 'Col. 5', '...');

%Print a horizontal line below the column headings.
fprintf('%s\n','---------------------------------------------------------------------------------------------------------');
fprintf(OutputFile, '%s\n','---------------------------------------------------------------------------------------------------------');

%Because Matlab arrays can't start at index 0, R is of size kmax+1 by kmax+1.
%R_0,0 (from the Romberg Integration handout) is stored in R(1,1).
%R_1,0 (from the Romberg Integration handout) is stored in R(2,1).
%R_1,1 (from the Romberg Integration handout) is stored in R(2,2).
%However, as mentioned in a note at the start of this code, you may want to rethink whether
%you really even want a square array in the first place.

R = zeros(1,kmax+1);

%Before entering the main loop, evaluate R(1,1) and print it.  Note that R(1,1) --- which is really
%R_0,0 in the Romberg table on the Romberg Integration handout --- is the only application of
%(Composite) Trapezoidal Rule that's NOT done via Eq. (5) of the handout.
R(k+1) = 0.5*h*(f(a) + f(b));
fprintf('  %2d %11.7f\n', k, R(1));
fprintf(OutputFile, '  %2d %11.7f\n', k, R(1));

%Main loop
while (check >= epsilon) && (k < kmax)
    
    Rold = R(1);
    Rcheck = R(k+1);
    
    %Apply Composite Trapezoidal Rule once, producing R_k+1,0
    sum = 0;
    for j = 1:2^k
        sum = sum + f(a + ((2*j)-1)*(h/(2^(k+1))));
    end
    R(1) = (0.5*Rold) + (h/(2^(k+1)))*sum;
    
    %compute the rest of the entries in the k+1st row of the table.
    for j = 1:k+1
        Rsave = R(j+1);
        R(j+1)= C(j)*(R(j)- Rold) + Rold;
        Rold = Rsave;
    end
    
    check = abs(R(k+1) - Rcheck);
    %Now print one row of the results table.
    fprintf('  %2d %11.7f', k+1, R(1));
    fprintf(OutputFile, '  %2d %11.7f', k+1, R(1));
    for j=2:k+2
        fprintf(' %11.7f', R(j));
        fprintf(OutputFile, ' %11.7f', R(j));
    end
    fprintf('\n');
    fprintf(OutputFile, '\n');
    
    %Increment the iteration counter.
    k = k + 1;
end

%Print another horizontal line.
fprintf('%s\n','---------------------------------------------------------------------------------------------------------');
fprintf(OutputFile, '%s\n','---------------------------------------------------------------------------------------------------------');

%Print a conclusion statement.
if check >= epsilon
    fprintf(' The method did not converge in %2d iterations.\n\n', kmax);
    fprintf(OutputFile, ' The method did not converge in %2d iterations.\n\n', kmax);
else
    fprintf(' The method converged to %12.14f after %2d iterations.\n\n', R(k+1), k);
    fprintf(OutputFile, ' The method converged to %12.14f after %2d iterations.\n\n', R(k+1), k);
end

Error = Actual - R(k+1);

F = 2;
for i = 1:k
    F = F + i;
end

Nf = F/(-log10(abs(Error)));

fprintf(' N_f = %12.7f .\n\n', Nf);
fprintf(OutputFile, ' N_f = %12.7f .\n\n', Nf);

%Close the output file.
fclose(OutputFile);

end

function C = GetCvalues(kmax)
C = zeros(1,kmax);

for j = 1:kmax
    C(j) = (4^j)/((4^j)-1);
end

%Compute all of the C values using Eq. (16) from the Romberg Integration handout.

return
end

% function answer = f(x)
% 
% answer = (x^3)*log(x);
% 
% return
% end

% function answer = f(x)
% 
% answer = (x^2*cos(3*x));
% 
% return
% end

% function answer = f(x)
% 
% answer = (sin(2*x)^2)/x;
% 
% return
% end