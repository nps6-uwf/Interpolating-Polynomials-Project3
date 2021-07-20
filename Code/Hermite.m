function [H a b]=Hermite(u,x,y,dy)  % Hermite interpolation (Lagrange)
                        % u: discrete data points;
                        % vector x: [x_1,...,x_n]
                        % vector y: [y_1,...,y_n]
                        % vector dy: [y'_1,...,y'_n]
  n=length(x);          % number of interpolating points
  k=length(u);          % number of discrete data points
  li=ones(n,k);         % Lagrange basis polynomials
  a=zeros(n,k);         % basis polynomials alpha(x)
  b=zeros(n,k);         % basis polynomials beta(x)   
  H=zeros(1,k);         % Hermie interpolation polynomial H(x)
  for i=1:n      
      dl=0;             % derivative of Lagrange basis
      for j=1:n   
          if j~=i
              dl=dl+1/(x(i)-x(j));
              li(i,:)=li(i,:).*(u-x(j))/(x(i)-x(j));
          end
      end
      l2=li(i,:).^2;
      b(i,:)=(u-x(i)).*l2;           % basis polynomial alpha(x)
      a(i,:)=(1-2*(u-x(i))*dl).*l2;  % basis polynomial beta(x)
      H=H+a(i,:)*y(i)+b(i,:)*dy(i);  % Hermite polynomial H(x)
  end  
end
