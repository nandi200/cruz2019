%%implementación del metodo númerico RK4, para la resolución de 
%% una ecuación diferencial

%%condiciones iniciales

function r=resolv

format long

x0 = 0; y = 1; x = 2; h = 0.2;

r=rungeKutta(x0,y,x,h);

