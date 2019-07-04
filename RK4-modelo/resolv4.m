%%implementación del metodo númerico RK4, para la resolución de 
%%sistema de dos ecuaciones diferenciales
%%utilizando una matriz para las ecuaciones
%%análisis de sensibilidad

%%%%%1. tomando las ctes de la literatura 
%%%%%2. Promoner rangos para los parametros a medir (-1% 1%)
%%%%%     a=-0.1;         ->[-0.101 -0.009]
%%%%%     b=0.1;          ->[0.009 0.101]
%%%%%     c=0.2;          ->[0.198 0.202]
%%%%%3. Evaluar el modelo con los parametros para cada rango se tomaran 50
%%%%%   muestras para cada rango
%%%%%En cada iteracion se calcula un error promedio con los valores (entre el
%%%%%valor de referencia y la simulacion con los parametros) -> dy/dyref NRSM


function r=resolv4
format long

    %%condiciones del RungeKutta
    tf= 19.5;%%%20 iteraciones iniciando en 0 hasta 19.5 
             %%%(20-0.5) cada iteracion vale 0.5
    h=0.5; % tamaño de paso
    
    %Archivo de salida
    fi = fopen('sensitivity_EDO2.txt','w+');
    
    %%condiciones iniciales
    c(1)=0;  %% tiempo inicial
    c(2)=5000;  %%x0 inicial
    c(3)=0; %%y0 inicial

    m=2; %%número de ecuaciones
    
    
    %%arreglos para almacenar los resultados
    Y=[];
    fh=[];
    fe=[];
   
    
    
    %%opcion para buscar sensibilidad de a,b,c
    op=1;
    switch op
        case 1
            refX=-0.1;
        case 2
            refX=0.1;
        case 3
            refX=0.2;
    end
            
    %%funcion de referencia, parametros normales  
    resulA=rungeKutta4(tf,c,h,m,op,refX),         
    
    rang=0.01;%%%%rango en que se evalua la sensibilidad +1%,-1%
    tamM=10;%%%tamaño del vector de muestra    
    rangI=refX-(rang)*(refX);
    rangF=refX+(rang)*(refX);

    %%numero de iteraciones por opcion
    X=linspace(rangI,rangF,tamM);
       
for i=1:length(X)
     %%condiciones del RungeKutta
    tf= 19.5;%%%20 iteraciones
    h=0.5; % tamaño de paso
     
    %%condiciones iniciales
    c(1)=0;  %% tiempo inicial
    c(2)=5000;  %%x0 inicial
    c(3)=0; %%y0 inicial

    m=2; %%número de ecuaciones
    
    %arreglos para almacenar los resultados
    Y=[];
    resulB=[];
    
    resulB=rungeKutta4(tf,c,h,m,op,X(i)),
     
     
    fh=[fh,resulB]; %%vector con las soluciones de la primera ecu 
    fe=[fe,X(i)];   %%vector con el parametro variado
end

op,

%%error cuadratico medio
ECM=[];
tam=size(fh);
for j=1:tam(2) %%columna
    suma=0;
    for i=1:tam(1) %fila
          suma=suma+(fh(i,j)-resulA(i))^2
    end
    ECM(j)=sqrt((suma)/tam(1));
end

sumaECM =sum(ECM);
sizeECM=size(ECM);
sensi1=sumaECM/sizeECM(2);

 sumaFe=sum(fe);
 sizeFe=size(fe);
 promFe=sumaFe/sizeFe(2);
 
 sensi2=abs((promFe-refX)/refX);
 if not(isreal(sensi2))|| isnan(sensi2)|| isinf(sensi2) || sensi2==0 || sensi2<1*10^-15 %%cuando son iguales resulA y resulB
     sensi2=1,
 end
 
 sensi=sensi1/sensi2,

fclose(fi);
end