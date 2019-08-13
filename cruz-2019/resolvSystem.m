%%implementaci�n del metodo n�merico RK4, para la resoluci�n de 
%%sistema de EDO
%%utilizando una matriz para las ecuaciones
%%an�lisis de sensibilidad

%%observaciones
%%simulaci�n de 0 a 365

%%erandi castillo



function r=resolvSystem
format long
    %%condiciones del RungeKutta
    tf= 19.5;%%%20 iteraciones iniciando en 0 hasta 19.5 
             %%%(20-0.5) cada iteracion vale 0.5
    h=0.5; % tama�o de paso
    
    %Archivo de salida
    fi = fopen('sensitivity.txt','w+'); 
    
    
    
 % % % % condiciones iniciales
        c(1)=0;    %tiempo inicial
        c(2)=832603;              %sh(0)
        c(3)=0.05*975641;         %ih(0)<--------falta esta condición inicial        
        c(4)=0.95*975641;         %sv
        c(5)=0.05*975641;         %iv(0)
        c(6)=866743               %ethc
        c(7)=866743               %enthc
        c(8)=317083               %lthc
        c(9)=317083               %lnthc
        c(10)=487820               %pthc
        c(11)=487820              %pnthc
    
     

    m=10; %%n�mero de ecuaciones
    
    
    %%arreglos para almacenar los resultados
    Y=[];
    fh=[];
    fe=[];
   
    
    %%%variables
    %%opcion para buscar sensibilidad en los parametros
    op=1;    
    
     switch op
        case 1
            refX=0.000037530493; %% muh=0.000037530493; -->fijo no variar
        case 2 
            refX=0.65;             %%a=0.65;   
        case 3
            refX=0.425;             %%bvh=0.425;  [0.1,0.75] day^-1
        case 4
            refX=0.425;             %%  bhv=0.425; [0.1,0.75] day^-1
        case 5
            refX=0.175;             %%ro=0.175; 1/7[1/10,1/4] day^-1
        case 6
            refX=0.07916666667;     %%muv [1/30,1/8] day^-1
        case 7 
            refX=3.5;             %%tthc [1,6] day^-1
        case 8 
            refX=3.5;             %%tnthc [1,6] day^-1
        case 9
           refX=50050;              %%emthc[10^3,10^6]
        case 10
            refX=50050;             %%emnthc [10^3,10^6]
        case 11
            refX=0.7;             %%ethc[0.7] day^-1
        case 12
            refX=0.7;             %%enthc[0.7] day^-1
        case 13
            refX=0.3;             %%muE[0.2,0.4]
        case 14
            refX=250250;             %%lmthc[5*10^2,5*10^5]
        case 15     
            refX=250250;             %%lmnthc[5*10^2,5*10^5]
        case 16
            refX=0.5;             %%lthc 0.5 day^-1
        case 17
            refX=0.5;             %%lnthc 0.5 day^-1
        case 18
            refX=0.3;             %%muL[0.2,0.4]
            
            %--->no variar muP,zi,fi
        case 19
            refX=0.115;             %%gv [0.08,0.15]
           %% gthc, gnthc son dependientes de zi,pi,t y fi
           %%no variar c0
        case 20
             refX=0.456;            %%kf 0.465, [0,0.571]
        case 21
             refX=0.014;            %%r 0.014,[0.010,0.018]
         
     end
    
    %%funcion de referencia, parametros normales  
    resulA=rungeKutta(tf,c,h,m,op,refX),         
    
    rang=0.1;%%%%rango en que se evalua la sensibilidad +1%,-1%
    tamM=10;%%%tama�o del vector de muestra    
    rangI=refX-(rang)*(refX);
    rangF=refX+(rang)*(refX);

    %%numero de iteraciones por opcion
    X=linspace(rangI,rangF,tamM);
for i=1:length(X)
     %%condiciones del RungeKutta
    tf= 19.5;%%%20 iteraciones
    h=0.5; % tama�o de paso
     
   
    % % % % condiciones iniciales
        c(1)=0;    %tiempo inicial
        c(2)=832603;              %sh(0)
        c(3)=0.05*975641;         %ih(0)<--------falta esta condición inicial        
        c(4)=0.95*975641;         %sv
        c(5)=0.05*975641;         %iv(0)
        c(6)=866743               %ethc
        c(7)=866743               %enthc
        c(8)=317083               %lthc
        c(9)=317083               %lnthc
        c(10)=487820               %pthc
        c(11)=487820              %pnthc
    

    m=10; %%n�mero de ecuaciones
    
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