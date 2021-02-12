function [A1,B1, A_1_error, B_1_error] = solve_line_fit(n, time,Temp, modif)
%Linus Schmitz 
%Final Version 10/25/2020
%Take in two vectors, one for time and one for tempurature 
t = time;
y = Temp;

% figure(n)
% hold on
% plot(t,y)

%Determine the length of either vector (They must be the same length)
N = length(t);

%Calculate the delta value(This is just the deominator of the first
%equaiton a)
delta = N * sum(t.^2) - sum(t)^2;
%Calculate first equaion (This is the same as A1 output and in the linear
%equation format is the b value which shifts the line in the vertical
%direction. 
a = (N * sum(t .* y) - sum(t) * sum (y))/ delta;
b = (1/N) * ( sum(y) - a * sum(t));

H = [ones(N,1),t];
sigma_y = sqrt((1/(N-2)) * sum((y - a -(b * t)).^2));
%The sigma_y value is far too large for it to be reasonable to similar to
%the coding challenge 5 we rewrote a hardcoded value for it instead. 
%sigma_y = 1;

W = diag(ones(1, N));
W =  W *(1/(sigma_y^2));

P = (inv(H' * W * H));
x_hat = (P * H' * W * y);

A1 = x_hat(1) + modif;
B1 = x_hat(2);

%B_1_error is the error in the preidcted value. The A_1_error is an
%arbitrary value for the error in the initial value which is not reasonable

A_1_error = sqrt(P(1));
%A_1_error = sigma_y * sqrt(sum(t.^2)/delta1);
B_1_error = sqrt(P(4));
%B_1_error = sigma_y * sqrt(N/delta1);


%Plot the new line to see if it is the same. 
% G1 = A1 + B1 .*t;
% plot(t,G1)
% legend('temps', 'G1' )
% hold off


end