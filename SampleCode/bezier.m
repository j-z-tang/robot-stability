function sum = bezier(a,s)
% Compute bezier polynomial output and derivative...
% Compute the Bezier polynomial
sum = 0;
order = length(a) - 1;
% compute bezier
for j = 1:(order+1)
    k = j-1;
    sum = sum + a(j)*(s^k)*((1-s)^(order-k))*factorial(order)/(factorial(k)*factorial(order-k));
end
end

