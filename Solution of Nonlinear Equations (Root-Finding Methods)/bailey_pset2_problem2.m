function bailey_pset2_problem2

%Clear the command window.
clc;
%Clear all previous variables.
clear all;
%Close all previously opened figures or images.
close all;

%Open an output file.
OutputFile = fopen('bailey_pset2_problem2.txt','w');

%Print title to the screen and to the output file.
fprintf('\n OUTPUT FROM bailey_pset2_problem2.m \n\n');
fprintf(OutputFile, '\n OUTPUT FROM bailey_pset2_problem2.m \n\n');

%Initialize endpoints of interval.
T = input(' Please input the temperature, T. \n');   %Left endpoint of interval
p = input(' Please input the pressure, p. \n');  %Right endpoint of interval
R = 0.08205;


v_0 = (R*T)/p;
v_1 = (R*(T))/p +.01;
v_2 = (R*T)/(p)+.03;

C_1 = 0.05587;
C_2 = 2.2769;
C_3 = 128300;
C_4 = -0.01587;
C_5 = 0.01855;

function answer = f(v)
B = (C_1*R*T) - C_2 - (C_3*R)/(T^2);
g = -(C_1*C_4*R*T) + (C_2*C_5) - (C_1*C_3*R)/(T^2);
D = (C_1*C_3*C_4*R)/(T^2);

answer = ((R*T)/v)+(B/(v^2))+(g/(v^3))+(D/(v^4))-p ;

return
end

%Initialize some other variables.
n = 0;                %Iteration counter
nmax = 50;            %Maximum allowed number of iterations
epsilon = 0.00001;    %Tolerance


fv_0 = f(v_0);
fv_1 = f(v_1);
fv_2 = f(v_2);


%Print information about the method and the problem to the screen and to the output file.
fprintf(' Solving a nonlinear equation using Muller’s Method \n');
fprintf(' with initial input T = %9.6f and p = %9.6f \n\n', T, p);
fprintf(' with initial guesses v_2 = %9.6f, v_1 = %9.6f and v_0= %9.6f \n\n', v_2, v_1, v_0);
fprintf(' with initial function vaules fv_2 = %9.6f, fv_1 = %9.6f and fv_0= %9.6f \n\n', fv_2, fv_1, fv_0);
fprintf(OutputFile, ' Solving a nonlinear equation using Muller’s Method \n');
fprintf(OutputFile, ' with initial input T = %9.6f and p = %9.6f \n\n', T, p);
fprintf(OutputFile, ' with initial guesses v_2 = %9.6f, v_1 = %9.6f and v_0= %9.6f \n\n', v_2, v_1, v_0);
fprintf(OutputFile, ' with initial function vaules fv_2 = %9.6f, fv_1 = %9.6f and fv_0= %9.6f \n\n', fv_2, fv_1, fv_0);

%Print the column headings for the results table.
fprintf('%10s%22s%28s\n', 'Iteration', 'x', 'f(x)');
fprintf(OutputFile, '%10s%22s%28s\n', 'Iteration', 'x', 'f(x)');

%Print a horizontal line below the column headings.
fprintf('%s\n','------------------------------------------------------------------');
fprintf(OutputFile, '%s\n','------------------------------------------------------------------');


check = 0;
%Main loop
while (check == 0) && (n < nmax)
    d = (v_2 - v_0)*(v_1 - v_0)*(v_2 - v_1);
    a = ((-(v_2 - v_0)*(fv_1 - fv_0)) + ((v_1 - v_0)*(fv_2 - fv_0)))/d;
    b =(((v_2 - v_0)^2*(fv_1 - fv_0)) - ((v_1 - v_0)^2*(fv_2 - fv_0)))/d;
    c = fv_0;
    if n+1 > nmax
        fprintf('The nmax has been reached \n')
        break
    elseif ((b^2)-4*a*c) < 0
        fprintf('A complex number has occured \n')
        break
    end
    if b >= 0
        sgnb = 1;
    else
        sgnb = -1;
    end
    vnew = v_0-((2*c)/(b+sgnb*sqrt((b^2)-4*a*c)));
    check = (abs(vnew - v_0)/abs(vnew))<epsilon;
    if check == 1
        break
    end
    fvnew = f(vnew);
    %Now we print the iteration number n+1, the value of xnew, and the value of check
    fprintf('    %2d                  %12.6f               %12.6f\n', n+1, vnew, fvnew);
    fprintf(OutputFile, '    %2d                  %12.6f               %12.6f\n', n+1, vnew, fvnew);
    fv_2 = fv_1;
    fv_1 = fv_0;
    fv_0 = fvnew;
    v_2 = v_1;
    v_1 = v_0;
    v_0 = vnew;
    n = n+1;
%     fprintf(' %12.6f\n',check);
end

%Print another horizontal line.
fprintf('%s\n','------------------------------------------------------------------');
fprintf(OutputFile, '%s\n','------------------------------------------------------------------');

%Print a conclusion statement.
if check == 1
    fprintf(' The method converged to %12.6f after %2d iterations.\n\n', vnew, n);
    fprintf(OutputFile, ' The method converged to %12.6f after %2d iterations.\n\n', vnew, n);
elseif (b^2-4*a*c) < 0
    fprintf(' The method stopped after %2d iterations because a complex number was going to be computed. \n\n', n);
    fprintf(OutputFile, 'The method stopped after %2d iterations because a complex number was going to be computed. \n\n', n);
else
    fprintf(' The method did not converge in %2d iterations.\n\n', n);
    fprintf(OutputFile, ' The method did not converge in %2d iterations.\n\n', n);
end

%Close the output file.
fclose(OutputFile);
end