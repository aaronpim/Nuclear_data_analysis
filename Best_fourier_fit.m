function [p,y_app] = Best_fourier_fit(x,y)
N = length(x);
p = zeros(N+1,2);

sine = @(m) sin(2*pi*m*(x-min(x))/range(x));
cosine = @(m) sin(2*pi*m*(x-min(x))/range(x));

coef_sin = @(m) (2/range(x))*trapz(y.*sine(m),x-min(x));
coef_cos = @(m) (2/range(x))*trapz(y.*cosine(m),x-min(x));

p(1,1) = 0.5*coef_cos(0);
y_app = p(1,1);
for n = 1:N
    p(n+1,1) = coef_cos(n);
    p(n+1,2) = coef_sin(n);
    y_app = y_app + p(n+1,1)*cosine(n) + p(n+1,2)*sine(n);
end