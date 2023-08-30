function [p,y_app] = Best_fourier_fit(x,y,alpha)
xs = min(x):min(diff(x)):max(x);
s = spline(x,y,xs);
B = fft(s);
[~,I] = maxk(abs(B),round(numel(x)*alpha)); 
ffts = zeros(length(B),1); ffts(I) = B(I); 
p = sparse(ffts);
y_app = spline(xs,real(ifft(ffts)),x);
end