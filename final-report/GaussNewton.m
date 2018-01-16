function [unknowns,steps,S] = GaussNewton()
format long
tolerance = 1e-8;
maxstep = 100;
p = rand(30,1);
q = rand(30,1);
d = rand(30,1);
a = [0.5;0.5];
m = length(p);
n = length(a);
aold = a;
for k = 1:maxstep
    S = 0;
    for i = 1:m
        for j = 1:n
            J(i,j) = df(p(i),q(i),a(1,1),a(2,1),j);
            JT(j,i) = J(i,j);
        end
    end
    Jz = -JT*J;
    for i = 1:m
        r(i,1) = d(i) - sqrt((a(1,1) - p(i))^2 + (a(2,1) - q(i))^2);
        S = S + r(i,1)^2;
    end
    k
    S
    g = Jz\JT;
    a = aold - g*r;
    unknowns = a;
    error(k) = a(1,1) - aold(1,1);
    if (abs(error(k)) <= tolerance)
        break;
    end
    aold = a;
end
steps = k;
hold all
plot (p,q,'r*')
plot(a(1,1),a(2,1),'b*')
title('Gauss-Newton Approximation of randomly placed shared bikes')
xlabel('X')
ylabel('Y')
legend('Data Points','Gauss-Newton Approximation of the user')
hold off
errorration(3:k) = error(2:k-1)./error(3:k);
%errorration
end

function value = df(p,q,a1,a2,index)
switch index
    case 1
        value = (2*a1 - 2*p)*0.5*((a1-p)^2 + (a2-q)^2)^(-0.5);
    case 2
        value = (2*a2 - 2*q)*0.5*((a1-p)^2 + (a2-q)^2)^(-0.5);
end
end
