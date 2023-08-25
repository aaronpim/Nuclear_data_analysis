A = importdata('ENDF_U_235_N_TOT_SIG.txt');
alpha = 0.001;
% This is data extracted from the followin web address
% https://www-nds.iaea.org/exfor/servlet/X4sShowData?db=x4&op=get_plotdata&req=-1&ii=5503&File=E4R18464_e4.zvd.dat.txt

x = log10(A(:,1)); y = log10(A(:,2));
[x,I] = unique(x); y = y(I);
% Convert the data into x and y coordinates, and eliminate any duplicate
% data.

dydx = diff(y)./diff(x);
R = abs(dydx)./max(abs(dydx));
% Identify areas of highly oscillatory behaviour.

subplot(2,1,1); plot(x,y); xlabel('position','Interpreter','latex'); ylabel('$\sigma$','Interpreter','latex');
subplot(2,1,2); plot(x(1:end-1),R<alpha); xlabel('position','Interpreter','latex'); ylabel('$\left|\frac{\partial \sigma}{\partial x}\right|< \alpha \times \left\|\frac{\partial \sigma}{\partial x}\right\|_{\infty}$','Interpreter','latex');
% Plot both side by side for ease of idenification.

x_poly_1 = x(1:find((R<alpha)==0,1)-1); y_1 = y(1:find((R<alpha)==0,1)-1); 
x_poly_2 = x(find((R<alpha)==0,1,"last")+1:length(x)); y_2 = y(find((R<alpha)==0,1,"last")+1:length(x));
x_ns = x(find((R<alpha)==0,1):find((R<alpha)==0,1,"last")); y_ns = y(find((R<alpha)==0,1):find((R<alpha)==0,1,"last"));
% Select the ranges of x-values for which a polynomial approximation might work.

%% Polynomial Approximation in the smooth regions
p1 = Best_poly_fit(x_poly_1,y_1); f1 = polyval(p1,x_poly_1);
p2 = Best_poly_fit(x_poly_2,y_2); f2 = polyval(p2,x_poly_2);
figure;
subplot(2,2,1); plot(x_poly_1,y_1,x_poly_1,f1);
subplot(2,2,2); plot(x_poly_2,y_2,x_poly_2,f2);
subplot(2,2,3); plot(x_poly_1,y_1 - f1);
subplot(2,2,4); plot(x_poly_2,y_2 - f2);
%% Fourier expansion in the non-smooth regions
p_ns = Best_poly_fit(x_ns,y_ns); f_ns = polyval(p_ns,x_ns); 
figure; subplot(2,1,1); plot(x_ns,y_ns,x_ns,f_ns);
g_ns = y_ns - f_ns;
subplot(2,1,2); plot(x_ns,g_ns);
[p,~] = Best_fourier_fit(x_ns,g_ns);
figure; histogram(p(:)); set(gca, 'Yscale', 'log');
xlabel('Fourier Coeficient'); ylabel('Frequency');