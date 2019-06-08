%%implementación del metodo númerico RK4, para la resolución de 
%% sistema de dos ecuaciones diferenciales

%%condiciones iniciales

function r=resolv2

format long


    tf= 19.5;%%%20 iteraciones iniciando en 0 hasta 19.5 
             %%%(20-0.5) cada iteracion vale 0.5
    h=0.5; % tamaño de paso
    c(1)=0;  %% tiempo inicial
    c(2)=5000;  %%x0 inicial
    c(3)=0; %%y0 inicial

    rungeKutta2(tf,c,h);

end
