%% Lagrange Interpolation. single vs double
clear all; clc;
f = @(x) sin(x); n=45; a = 0; b = 2*pi;
xts = single(linspace(a,b,1000));% test points
yts = f(xts); 
xtd = linspace(a,b,1000);% test points
ytd = f(xtd); 
for i=[6,11,16,26,36]
    xd = linspace(a,b,i);%interpolation points
    xs = single(linspace(a,b,i));%interpolation points
    yd = f(xd);% y values
    ys = f(xs);% y values
    yld = lagrangepoly(xd,yd,xtd);% Lagrange polynomial
    yls = lagrangepoly(xs,ys,xts);% Lagrange polynomial
    figure(1)
    semilogy(xts,abs(yts-yls),'-b',xtd,abs(ytd-yld),'-r')
    xlabel('x');    ylabel('absolute error');
    title(sprintf('Error of Lagrange Interpolation, #%i',i-1))
    legend('Lagrange-single','Lagrange-double')
    figure(2)
    plot(xts,yls,'-b',xtd,yld,'-r',xd,yd,'ok',xtd,ytd,'-k')
    xlabel('x');    ylabel('y');
    title(sprintf('Lagrante Interpolation, #%i',i-1))
    legend('Lagrange-single','Lagrange-double'); pause
end 
%% Direct Interpolation
clear all; clc;
f = @(x) sin(x); nt = 1000;
a= 0; b = 2*pi;
xt = linspace(a,b,nt);% test points
yt = f(xt);
n=[5:5:50]; j=0;
for i=n
    j=j+1;
    x = linspace(a,b,i);%interpolation points
    y = [f(x)]';% y values
    V = vander(x);
    cond(V)
    % Gaussain elimination
    pg=V\y;
    % inverse of matrix
    pinv = inv(V)*y;
    % Eigen-Decomposition
    [S,D] = eig(V);
    pe = S*((inv(S)*y)./diag(D));
    % SVD
    [U,S,X] = svd(V);
    ps = X*((U'*y)./diag(S));
    % QR
    [Q,R] = qr(V);
    pq =  R\(Q'*y);
    %
    yg = polyval(pg,xt);
    yi = polyval(pinv,xt);
    ye = polyval(pe,xt);
    ys = polyval(ps,xt);
    yq = polyval(pq,xt);
    %
    MaxErrg(j) = norm(yt-yg);
    MaxErri(j) = norm(yt-yi);
    MaxErre(j) = norm(yt-ye);
    MaxErrs(j) = norm(yt-ys);
    MaxErrq(j) = norm(yt-yq);
    %
    figure(1)
    semilogy(xt,abs(yt-yg),'-b',xt,abs(yt-yi),'-r',xt,abs(yt-ye),'--b',...
        xt,abs(yt-ys),'--r',xt,abs(yt-yq),'-k')
    xlabel('x');    ylabel('absolute error');
    title(sprintf('Error of Direct Method, #%i',i-1))
    legend('GE','Inv','Eig','SVD','QR')
    figure(2)
    plot(xt,yg,'-b',xt,yi,'-r',xt,ye,'--b',xt,ys,'--r',xt,yq,'-k',x,y,'ok',xt,yt,'-k')
    xlabel('x');    ylabel('y');
    title(sprintf('Direct Method, #%i',i-1))
    legend('GE','Inv','Eig','SVD','QR');
    pause(0.1)
end
%
figure(3)
semilogy(n,MaxErrg,'-b',n,MaxErri,'-r',n,MaxErre,'--b',n,MaxErrs,'--r',n,MaxErrq,'-k')
xlabel('x');    ylabel('absolute error');
title('Error of Direct Method')
legend('GE','Inv','Eig','SVD','QR')%%
% you created a figure and it is "current". % the following, you could have guessed
% set(gcf,'color','none');
% set(gca,'color','none');
% but this following little detail took me ages to figure out
% set(gcf,'InvertHardCopy','off');
%% Direct Interpolation
clear all; clc;
f = @(x) sin(x); nt = 1000;
a= 0; b = 2*pi;
xt = linspace(a,b,nt);% test points
yt = f(xt);
n=[5:5:100]; j=0;
for i=n
    j=j+1;
    x = linspace(a,b,i);%interpolation points
    y = [f(x)]';% y values
    V = vander(x);
    tic
    % Gaussain elimination
    pg=V\y;
    tg(j) = toc; tic
    % inverse of matrix
    pinvV = inv(V)*y;
    ti(j) = toc; tic
    % Eigen-Decomposition
    [S,D] = eig(V);
    pe = S*((inv(S)*y)./diag(D));
    te(j) = toc; tic
    % SVD
    [U,S,X] = svd(V);
    ps = X*((U'*y)./diag(S));
    ts(j) = toc; tic
    % QR
    [Q,R] = qr(V);
    pq =  R\(Q'*y);
    tq(j) = toc;
end
%
figure(1)
semilogy(n,tg,'-b',n,ti,'-r',n,te,'--b',n,ts,'--r',n,tq,'-k')
xlabel('Number of Interpolation Points');    ylabel('Time');
title(sprintf('CPU Time'))
legend('GE','Inv','Eig','SVD','QR')
