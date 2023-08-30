A = importdata('ENDF_U_235_N_TOT_SIG.txt');
alpha = 0.0025;
% This is data extracted from the followin web address
% https://www-nds.iaea.org/exfor/servlet/X4sShowData?db=x4&op=get_plotdata&req=-1&ii=5503&File=E4R18464_e4.zvd.dat.txt

x = log10(A(:,1)); y = log10(A(:,2));
[x,I] = unique(x); y = y(I); clearvars I A;
% Convert the data into x and y coordinates, and eliminate any duplicate
% data.

dydx = diff(y)./diff(x);
R = abs(dydx)./max(abs(dydx));
x1 = x(1:find((R<alpha)==0,1)-1); y1 = y(1:find((R<alpha)==0,1)-1); 
x2 = x(find((R<alpha)==0,1,"last")+1:length(x)); y2 = y(find((R<alpha)==0,1,"last")+1:length(x));
x_ns = x(find((R<alpha)==0,1):find((R<alpha)==0,1,"last")); y_ns = y(find((R<alpha)==0,1):find((R<alpha)==0,1,"last"));

%% Polynomial Approximation in the smooth regions
p1 = Best_poly_fit(x1,y1); f1 = polyval(p1,x1);
p2 = Best_poly_fit(x2,y2); f2 = polyval(p2,x2);

%% Approximate the remainder by a Fourier series
[~,s1] = Best_fourier_fit(x1,y1-f1,0.2);
[~,s2] = Best_fourier_fit(x2,y2-f2,0.2);
s1 = s1(2:end-1); s2 = s2(2:end-1); f1 = f1(2:end-1); f2 = f2(2:end-1); 
x1 = x1(2:end-1); x2 = x2(2:end-1); y1 = y1(2:end-1); y2 = y2(2:end-1); 
%% Approximate the non-smooth section
[~,s_ns] = Best_fourier_fit(x_ns,y_ns,0.2);
s_ns = s_ns(2:end-1); x_ns = x_ns(2:end-1); y_ns = y_ns(2:end-1);
%% Plot the function
figure; plot(x1,f1+s1,'b-',x2,f2+s2,'b-',x_ns,s_ns,'b-',x,y,'r--'); xlabel('Position'); ylabel('Cross-section')
figure; plot(x1,(f1+s1-y1)./abs(y1),'k-',x2,(f2+s2-y2)./abs(y2),'k-',x_ns,(s_ns-y_ns)./abs(y_ns),'b-'); xlabel('Position'); ylabel('Relative Error')