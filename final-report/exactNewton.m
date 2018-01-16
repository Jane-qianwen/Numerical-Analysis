syms u v
format long
tolerance = 1e-8;
initial_estimate = [0.5,0.5];
maxstep = 100;
p = rand(30,1);
q = rand(30,1);
d = rand(30,1);
m = length(p);
n = length(initial_estimate);
F1 = d(1) - sqrt((u - p(1))^2 + (v - q(1))^2);
F2 = d(2) - sqrt((u - p(2))^2 + (v - q(2))^2);
F3 = d(3) - sqrt((u - p(3))^2 + (v - q(3))^2);
F4 = d(4) - sqrt((u - p(4))^2 + (v - q(4))^2);
F5 = d(5) - sqrt((u - p(5))^2 + (v - q(5))^2);
F6 = d(6) - sqrt((u - p(6))^2 + (v - q(6))^2);
F7 = d(7) - sqrt((u - p(7))^2 + (v - q(7))^2);
F8 = d(8) - sqrt((u - p(8))^2 + (v - q(8))^2);
F9 = d(9) - sqrt((u - p(9))^2 + (v - q(9))^2);
F10 = d(10) - sqrt((u - p(10))^2 + (v - q(10))^2);
F11 = d(11) - sqrt((u - p(11))^2 + (v - q(11))^2);
F12 = d(12) - sqrt((u - p(12))^2 + (v - q(12))^2);
F13 = d(13) - sqrt((u - p(13))^2 + (v - q(13))^2);
F14 = d(14) - sqrt((u - p(14))^2 + (v - q(14))^2);
F15 = d(15) - sqrt((u - p(15))^2 + (v - q(15))^2);
F16 = d(16) - sqrt((u - p(16))^2 + (v - q(16))^2);
F17 = d(17) - sqrt((u - p(17))^2 + (v - q(17))^2);
F18 = d(18) - sqrt((u - p(18))^2 + (v - q(18))^2);
F19 = d(19) - sqrt((u - p(19))^2 + (v - q(19))^2);
F20 = d(20) - sqrt((u - p(20))^2 + (v - q(20))^2);
F21 = d(21) - sqrt((u - p(21))^2 + (v - q(21))^2);
F22 = d(22) - sqrt((u - p(22))^2 + (v - q(22))^2);
F23 = d(23) - sqrt((u - p(23))^2 + (v - q(23))^2);
F24 = d(24) - sqrt((u - p(24))^2 + (v - q(24))^2);
F25 = d(25) - sqrt((u - p(25))^2 + (v - q(25))^2);
F26 = d(26) - sqrt((u - p(26))^2 + (v - q(26))^2);
F27 = d(27) - sqrt((u - p(27))^2 + (v - q(27))^2);
F28 = d(28) - sqrt((u - p(28))^2 + (v - q(28))^2);
F29 = d(29) - sqrt((u - p(29))^2 + (v - q(29))^2);
F30 = d(30) - sqrt((u - p(30))^2 + (v - q(30))^2);

%a = size([u,v]);
%disp(a);
a = newton_n_dim(tolerance,initial_estimate,[u,v],[F1;F2;F3;F4;F5;F6;F7;F8;F9;F10;F11;F12;F13;F14;F15;F16;F17;F18;F19;F20;F21;F22;F23;F24;F25;F26;F27;F28;F29;F30]);
%a = newton_n_dim(tolerance,initial_estimate,[u,v],[F1;F2]);
solution = double(a);
disp(solution);

hold all
plot (p,q,'r*')
plot(solution(1,1),solution(1,2),'b*')
title('exact Newton Approximation of randomly placed shared bikes')
xlabel('X')
ylabel('Y')
legend('Data Points','exact Newton Approximation of the user')
hold off


function [X] = newton_n_dim(tolerance_rss,initial_estimate,sym_variables,sym_equations)

H = jacobian(sym_equations,sym_variables);
X = initial_estimate;

n_equations = 0;
if length(sym_equations)==length(sym_variables)
    n_equations = 1;
end

stop = 0;
time = 0;
while ~stop
        F_X = subs(sym_equations,sym_variables,X); 
        F_prime_X = subs(H,sym_variables,X);
        if ~isnumeric(F_prime_X)
            F_prime_X = eval(F_prime_X);
        end
    if n_equations ==1
        d_X = (F_prime_X^-1)*F_X;
    else %overdetermined solution, use generalized inverse matrix
        twice_difference = (F_prime_X.'*F_prime_X);
        d_X = (twice_difference\F_prime_X.')*F_X;
    end
    residual = 0;
    m = length(sym_equations);
    for i = 1:m
        residual = residual+F_X(i);
    disp(double(residual));
    end
    %disp(d_X);
    Xold = X;
    %disp(Xold);
    X = X - d_X.';
    %disp(X);
    error = abs(X(1,1)-Xold(1,1)); %calculate the least-square error
    if (error < tolerance_rss)
        stop = 1;
    end
    time = time + 1
end
end
