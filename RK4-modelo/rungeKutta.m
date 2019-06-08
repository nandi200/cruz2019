%%encuentra el valor de Y para cada valor de X
%%usando un tamaño de paso h y con valoles
%%iniciales y0 y x0

function y = rungeKutta(x0,y0,x,h)

%cuenta el número de iteraciones usando un tamaño de paso h
n=(x-x0)/h;

%itera el número de iteraciones calculadas
y=y0;

for i=1:n
    %%aplica el Runge Kutta para
    %%encontrar el siguiente valor de y
    k1=h*dxdy(x0,y);
    k2 = h*dxdy(x0 + 0.5*h, y + 0.5*k1);
    k3 = h*dxdy(x0 + 0.5*h, y + 0.5*k2);
    k4 = h*dxdy(x0 + h, y + k3);
    
    %%actualiza el valor de y
    y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
      
    %%actualiza el siguiente valor de x
    x0 = x0 + h;

end

end

