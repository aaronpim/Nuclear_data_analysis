function [x_transition,p1,four1,p2,four2,four_ns] = Nuclear_data_analysis_final(alpha, compress_smooth,compress_non_smooth,data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUTS
%
%   alpha                   a double-value scalar in the range 0 to 1 that 
%                           determines the tolerance for detecting the 
%                           transition between smooth regions of the data 
%                           and non-smooth.
%
%   compress_non_smooth     a double-value scalar in the range 0 to 1 that 
%                           determines the compression level of the fourier
%                           coefficients in the non-smooth region.
%
%   compress_smooth         a double-value scalar in the range 0 to 1 that 
%                           determines the compression level of the fourier
%                           coefficients in the smooth region.
%
%   data                    a string that is the file name of the data
%                           input. It is assumed that the data is two
%                           columns, seperated by spaces.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   OUTPUTS
%
%   x_transition            a vector of two values that determine the
%                           transition points between the smooth region and
%                           the non-smooth region.
%
%   p1                      a vector of polynomial coefficients that
%                           determine the describe the behaviour in the
%                           left smooth tail region of the data.
%
%   p2                      a vector of polynomial coefficients that
%                           determine the describe the behaviour in the
%                           right smooth tail region of the data.
%
%   four1                   a sparse vector that contains the compressed
%                           fourier coefficients of the left smooth tail 
%                           region of the data after the polynomial
%                           approximation has been subtracted.
%
%   four2                   a sparse vector that contains the compressed
%                           fourier coefficients of the right smooth tail 
%                           region of the data after the polynomial
%                           approximation has been subtracted.
%
%   four_ns                 a sparse vector that contains the compressed
%                           fourier coefficients of the non-smooth region
%                           of the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults
if nargin == 3
    data = 'ENDF_U_235_N_TOT_SIG.txt';
elseif nargin == 2
    data = 'ENDF_U_235_N_TOT_SIG.txt';
    compress_non_smooth = 0.2;
elseif nargin == 1
    data = 'ENDF_U_235_N_TOT_SIG.txt';
    compress_non_smooth = 0.2;
    compress_smooth = 0.2;
elseif nargin == 0
    data = 'ENDF_U_235_N_TOT_SIG.txt';
    compress_non_smooth = 0.2;
    compress_smooth = 0.2;
    alpha = 0.0025;
end

%% Importing and cleaning the data
A = importdata(data);
% This is data extracted from the followin web address
% https://www-nds.iaea.org/exfor/servlet/X4sShowData?db=x4&op=get_plotdata&req=-1&ii=5503&File=E4R18464_e4.zvd.dat.txt

x = log10(A(:,1)); y = log10(A(:,2));
[x,I] = unique(x); y = y(I); clearvars I A;
% Convert the data into x and y coordinates, and eliminate any duplicate
% data.

%% Determining transition of smooth to non-smooth.
dydx = diff(y)./diff(x);
R = abs(dydx)./max(abs(dydx));
% Calculate the derivative, when the derivative changes sharply then we
% consider that the "cut-off" point.

x1 = x(1:find((R<alpha)==0,1)-1); 
y1 = y(1:find((R<alpha)==0,1)-1); 
% The left-smooth tail region

x2 = x(find((R<alpha)==0,1,"last")+1:length(x)); 
y2 = y(find((R<alpha)==0,1,"last")+1:length(x));
% The right-smooth tail region

x_ns = x(find((R<alpha)==0,1):find((R<alpha)==0,1,"last")); 
y_ns = y(find((R<alpha)==0,1):find((R<alpha)==0,1,"last"));
% The non-smooth centre region

x_transition = [find((R<alpha)==0,1),find((R<alpha)==0,1,"last")];
% Note the transition points.

%% Polynomial Approximation in the smooth regions
p1 = Best_poly_fit(x1,y1); f1 = polyval(p1,x1);
p2 = Best_poly_fit(x2,y2); f2 = polyval(p2,x2);

%% Approximate the remainder by a Fourier series
[four1,s1] = Best_fourier_fit(x1,y1-f1,compress_smooth); four1 = sparse(four1);
[four2,s2] = Best_fourier_fit(x2,y2-f2,compress_smooth); four2 = sparse(four2);
s1 = s1(2:end-1); s2 = s2(2:end-1); f1 = f1(2:end-1); f2 = f2(2:end-1); 
x1 = x1(2:end-1); x2 = x2(2:end-1); y1 = y1(2:end-1); y2 = y2(2:end-1); 

%% Approximate the non-smooth section
[four_ns,s_ns] = Best_fourier_fit(x_ns,y_ns,compress_non_smooth);
s_ns = s_ns(2:end-1); x_ns = x_ns(2:end-1); y_ns = y_ns(2:end-1);

%% Plot the function
subplot(2,1,1); plot(x1,f1+s1,'b-',x2,f2+s2,'b-',x_ns,s_ns,'b-',x,y,'r--'); xlabel('Position'); ylabel('Cross-section')
subplot(2,1,2); plot(x1,(f1+s1-y1)./abs(y1),'k-',x2,(f2+s2-y2)./abs(y2),'k-',x_ns,(s_ns-y_ns)./abs(y_ns),'b-'); xlabel('Position'); ylabel('Relative Error')