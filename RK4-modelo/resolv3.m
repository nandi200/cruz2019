%%implementaci�n del metodo n�merico RK4, para la resoluci�n de 
%%sistema de dos ecuaciones diferenciales
%%utilizando una matriz para las ecuaciones

%%condiciones iniciales

function r=resolv3
format long


    tf= 19.5;%%%20 iteraciones iniciando en 0 hasta 19.5 
             %%%(20-0.5) cada iteracion vale 0.5
    h=0.5; % tama�o de paso
    c(1)=0;  %% tiempo inicial
    c(2)=5000;  %%x0 inicial
    c(3)=0; %%y0 inicial

    m=2; %%n�mero de ecuaciones
    Y=[];
    Y=rungeKutta3(tf,c,h,m),

end
