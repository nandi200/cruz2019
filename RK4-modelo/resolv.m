%%implementaci�n del metodo n�merico RK4, para la resoluci�n de 
%% una ecuaci�n diferencial

%%condiciones iniciales

function r=resolv

format long

x0 = 0; y = 1; x = 2; h = 0.2;

r=rungeKutta(x0,y,x,h);

