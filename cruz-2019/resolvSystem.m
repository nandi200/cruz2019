%%implementación del metodo númerico RK4, para la resolución de 
%%sistema de EDO
%%utilizando una matriz para las ecuaciones
%%análisis de sensibilidad

%%observaciones
%%simulación de 0 a 365



function r=resolvSystem
format long
    %%condiciones del RungeKutta
    tf= 19.5;%%%20 iteraciones iniciando en 0 hasta 19.5 
             %%%(20-0.5) cada iteracion vale 0.5
    h=0.5; % tamaño de paso
    
    %Archivo de salida
    fi = fopen('sensitivity.txt','w+');
    
     % % % % condiciones iniciales
        c(1)=;    %tiempo inicial
        c(2)=;    
        c(3)=;   
        c(4)=; 		
        c(5)=; 
        c(6)=;    
        c(7)=;    
        c(8)=;    
        c(9)=;    
        c(10)=;   
        c(11)=;  
     

    m=10; %%número de ecuaciones
    
    
    %%arreglos para almacenar los resultados
    Y=[];
    fh=[];
    fe=[];
   
    
    %%%variables
    %%opcion para buscar sensibilidad en los parametros
    op=;
   
     switch op
        case 1
            refX=1/(75*365)^1; %%muh day^2  -->fijo no variar
        case 2 
            refX=;             %%alpha[0.3,1] day^-1  
        case 3
            refX=;             %%bvh [0.1,0.75] day^-1
        case 4
            refX=;             %%bhv [0.1,0.75] day^-1
        case 5
            refX=;             %%ro 1/7[1/10,1/4] day^-1
        case 6
            refX=;             %%muv [1/30,1/8] day^-1
        case 7 
            refX=;             %%tthc [1,6] day^-1
        case 8 
            refX=;             %%tnthc [1,6] day^-1
        case 9
           refX=;              %%eMthc[10^3,10^6]
        case 10
            refX=;             %%eMtnthc [10^3,10^6]
        case 11
            refX=;             %%ethc[0.7] day^-1
        case 12
            refX=;             %%enthc[0.7] day^-1
        case 13
            refX=;             %%mue[0.2,0.4]
        case 14
            refX=;             %%lMthc[5*10^2,5*10^5]
        case 15     
            refX=;             %%lMthc[5*10^2,5*10^5]
        case 16
            refX=;             %%lathc 0.5 day^-1
        case 17
            refX=;             %%lanthc 0.5 day^-1
        case 18
            refX=;             %%mul[0.2,0.4]
        case 19
            refX=;             %%mup[0.4] --->no variar
         case 20
            refX=;             %%zeta 0.5
        case 21
            refX=;             %%fi pi/2
        case 22
            refX=;             %%gv [0.08,0.15]
        case 23
            refX=;             %%gthc =gv*(1+zeta cos(2pi/365)t+fi)
        case 24
            refX=;             %%gnthc=gv*(1+zeta cos(2pi/365)t+fi)
        case 25
            refX=;             %%C0=0.001-->fijo no variar
        case 26
             refX=;            %%kf 0.465, [0,0.571]
        case 27
             refX=;            %%r 0.014,[0.010,0.018]
         case 28
             refX=;           %%((kf*C0*exp^(rt))/(kf+C0*(exp(rt)-1))+(((0.571-kf)*C0*exp^(rt))/((0.571-kf)+((C0*(exp^(rt))-1))))))
         case 29
             refX=;           %%nh 832603
     end
    
    %%funcion de referencia, parametros normales  
    resulA=rungeKutta(tf,c,h,m,op,refX),         
    
    rang=0.1;%%%%rango en que se evalua la sensibilidad +1%,-1%
    tamM=10;%%%tamaño del vector de muestra    
    rangI=refX-(rang)*(refX);
    rangF=refX+(rang)*(refX);

    %%numero de iteraciones por opcion
    X=linspace(rangI,rangF,tamM);
for i=1:length(X)
     %%condiciones del RungeKutta
    tf= 19.5;%%%20 iteraciones
    h=0.5; % tamaño de paso
     
   
    % % % % condiciones iniciales
        c(1)=0;    %tiempo inicial
        c(2)=832603;              %sh(0)
        c(3)=0.95*975641;         %sv
        c(4)=0.05*975641;         %lv(0)
        c(5)=866743               %ethc
        c(6)=866743               %enthc
        c(7)=317083               %lthc
        c(8)=317083               %lnthc
        c(9)=487820               %pthc
        c(10)=487820              %pnthc
    

    m=10; %%número de ecuaciones
    
    %arreglos para almacenar los resultados
    Y=[];
    resulB=[];
    
    resulB=rungeKutta(tf,c,h,m,op,X(i)),
     
     
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