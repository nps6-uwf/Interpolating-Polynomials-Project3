%% Combination with Chebyshev nodes and equal spaced interpolation points.
clear all; clc;
func = 2;
switch func
    case 1
        f = @(x) sin(x);
        fp = @(x) cos(x);
        a = 0; b = 2*pi;
    case 2
        f = @(x) sin(5*x);
        fp = @(x) 5*cos(5*x);
        a = 0; b = 2*pi;
    case 3
        f = @(x) 1./(x.^2 + 1);
        fp = @(x) -(2.*x)./(x.^2 + 1).^2;
        a = -5; b = 5;
    otherwise
        disp('Undefine function selected.')
end
%
n=45;
xt = linspace(a,b,1000);% test points
yt = f(xt);
%
for i=2:n
    x = linspace(a,b,i);%interpolation points
    xc = (b-a)/2*cos((2*(1:i)-1)*pi/(2*i))+(b+a)/2; % Chebyshev interpolation points
    y = [f(x)]';% y values
    yc = [f(xc)]';% y values
    % Lagrange-equal spaced
    yl = lagrangepoly(x,y,xt);% Lagrange polynomial 
    % Lagrange-Chebyshev
    ylc = lagrangepoly(xc,yc,xt);% Lagrange polynomial
    % direct-equal spaced
    V = vander(x);
    p=V\y;% centers for Newton form:
    yd = polyval(p,xt); 
    % direct-Chebyshev
    V = vander(xc);
    p=V\yc;% centers for Newton form:
    ydc = polyval(p,xt); 
    % Hermite-equal spaced
    yp = fp(x);% y prime values
    [yh, ~, ~]=Hermite(xt,x,y,yp);
    % Hermite-Chebyshev
    yp = fp(xc);% y prime values
    [yhc, ~, ~]=Hermite(xt,xc,yc,yp);
    %
    MaxErr_le(i-1) = norm(yt-yl,inf);
    MaxErr_lc(i-1) = norm(yt-ylc,inf);
    MaxErr_de(i-1) = norm(yt-yd,inf);
    MaxErr_dc(i-1) = norm(yt-ydc,inf);
    MaxErr_he(i-1) = norm(yt-yh,inf);
    MaxErr_hc(i-1) = norm(yt-yhc,inf);  
    %
    figure(1)
    semilogy(xt,abs(yt-yl),'-b', xt,abs(yt-ylc),'--b',xt,abs(yt-yd),'-r',...
        xt,abs(yt-ydc),'--r',xt,abs(yt-yh),'-k',xt,abs(yt-yhc),'--k')
    xlabel('x');    ylabel('absolute error');
    title(sprintf('Error of Interpolations, #%i',i))
    legend('Lagrange-e','Lagrange-c','Direct-e','Direct-c','Hermite-e','Hermite-c');
    figure(2)
    plot(xt,yl,'-b',xt,ylc,'--b',xt,yd,'-r',xt,ydc,'--r',xt,yh,'-k',xt,yhc,'--k',xc,yc,'og',x,y,'ok',xt,yt,'-g')
    xlabel('x');    ylabel('y');
    title(sprintf('Comparison, #%i',i))
    legend('Lagrange-e','Lagrange-c','Direct-e','Direct-c','Hermite-e','Hermite-c');
    pause(0.1)
end
%
datasite=2:n;
figure(3)
semilogy(datasite,MaxErr_le,'-b',datasite,MaxErr_lc,'--b',...
    datasite,MaxErr_de,'-r',datasite,MaxErr_dc,'--r',...
    datasite,MaxErr_he,'-k',datasite,MaxErr_hc,'--k')
xlabel('Number of Interpolation Points, N');    ylabel('Maximum Error,\epsilon');
legend('Lagrange-e','Lagrange-c','Direct-e','Direct-c','Hermite-e','Hermite-c');
title(sprintf('Error of Interpolations, #%i',i))