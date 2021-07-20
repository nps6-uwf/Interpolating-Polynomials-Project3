% Project 3.
% Name: Nick Sebasco
% Code from Dr. Kuo Lecture 7/13/2021 - Graduate Numerical Analysis
% Version 2: + added my automated figure saving code that I wrote for 
% project 2.

clear all; clc;

% styles
bgColor = [0.8,0.8,0.8];
funcNames = {'Sin(x)','Sin(5x)','1/(x^2+1)'};

% Specify output file extension and output directory name.
figureFileExtension = ".fig"; % jpg, fig ...
figureDirectory = 'P3Figures';

% Create new directory to store figures.
if ~exist(figureDirectory, 'dir')
       mkdir(figureDirectory)
else
    % wipe all preexitsing files from the figureDirectory
    delete(join([figureDirectory,'\*']));
end

for func = [1,2,3]
    switch func
        case 1 % for f(x)
            f = @(x) sin(x); fp = @(x) cos(x);
            a = 0; b = 2* pi;
        case 2 % for g(x)
            f = @(x) sin(5*x); fp = @(x) 5*cos(5*x);
            a = 0; b = 2* pi;
        case 3 % for h(x)
            f = @(x) 1./(x.^2+1); fp = @(x) -(2*x)./(x.^2 + 1).^2;
            a = -5; b = 5;
        otherwise
            disp('Invalid func, please try again.')
    end

    delta = 0; % An arbitrary delta can be added to each i (optional).
    nt = 1000;
    xt = linspace(a,b,nt);
    yt = f(xt);

    for i=[5,10,15,25,35]
        i = i + delta;
        x = linspace(a,b,i);
        y = [f(x)]';
        xc = (b-a)/2*cos((2*(1:i)-1)*pi/(2*i))+(b+a)/2; % Chebyshev nodes
        yc = [f(xc)]';

        % (1)direct-equal
        V = vander(x);
        p = V\y;
        yd = polyval(p,xt);

        % (2)direct-chebyshev
        V = vander(xc);
        p = V\yc;
        ydc = polyval(p,xt);  

        % (3)lagrange-equal
        yl = lagrangepoly(x,y,xt);

        % (4)lagrange-chebyshev
        ylc = lagrangepoly(xc,yc,xt);

        % (5)hermite-equal
        yp = fp(x);
        [yh,~,~] = Hermite(xt,x,y,yp);

        % (6)hermite-chebyshev
        ypc = fp(xc);
        [yhc,~,~] = Hermite(xt,xc,yc,ypc);

        figure(1) % Approximation plots
        % set(gca,'Color',[0.8,0.8,0.8])  
        plot(xt,yd,'-b',xt,ydc,'--b',xt,yl,'-r',xt,ylc,'--r',xt,yh,'-g',xt,yhc,'--g',...
            x,y,'*b',xc,yc,'or',xt,yt,'k-')
        set(gca,'Color',bgColor) 
        xlabel('X')
        ylabel('Y')
        title(sprintf('%s, Interpolation, n = %d',funcNames{func},i))
        legend('Direct-E','Direct-C','Lagrange-E','Lagrange-C','Hermite-E','Hermite-C')
        axis tight;

        figure(2) % Error plots
        semilogy(xt,abs(yd-yt),'-b',xt,abs(ydc-yt),'--b',xt,abs(yl-yt),'-r',xt,abs(ylc-yt),...
            '--r',xt,abs(yh-yt),'-g',xt,abs(yhc-yt),'--g')
        set(gca,'Color',bgColor) 
        xlabel('X')
        ylabel('Absolute Error')
        title(sprintf('%s, Absolute Error, n = %d',funcNames{func},i))
        legend('Direct-E','Direct-C','Lagrange-E','Lagrange-C','Hermite-E','Hermite-C')
        axis tight;
        % pause
        
        % save figures to the output directory.
        saveas(figure(1),[pwd,sprintf('%s/function%d_%d%s',figureDirectory,func,i,figureFileExtension)]); 
        saveas(figure(2),[pwd,sprintf('%s/function%d_absoluteError_%d%s',figureDirectory,func,i,figureFileExtension)]);

    end
end

