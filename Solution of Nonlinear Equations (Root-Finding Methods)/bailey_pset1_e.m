function bailey_pset1_e
x=[-2:.05:8]
function answer = f(x)
answer = x.^2-4.*x.*(1 + cos(2.*x)) + 2;

return
end
plot(x,f(x))
yline(0,'-.k')
title('x^2 - 4x(1 + cos(2x)) + 2')
xlabel('x')
ylabel('f(x)')
end

