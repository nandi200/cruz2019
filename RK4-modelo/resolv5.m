%%implementación del metodo númerico RK4, para la resolución de 
%%sistema de diez ecuaciones diferenciales
%%utilizando una matriz para las ecuaciones
%%análisis de sensibilidad

%%%%%1. tomando las ctes de la literatura 
%%%%%2. Promoner rangos para los parametros a medir (-1% 1%)
% % % % muh     0.0000365296803652968	0.0000361643835616438	0.0000368949771689498	
% % % % alpha	1                           0.99                    1.01	
% % % % betavh	[0.1,0.75]                  0.4158                  0.4242              0.42
% % % % betahv	[0.1,0.75]                  0.4158                  0.4242              0.42
% % % % ro      0.1428                      0.141372                0.144228	
% % % % delta	0.001                       0.00099                 0.00101	
% % % % gammathc	0.08                    0.0792                  0.0808	
% % % % gammanthc	0.08                    0.0792                  0.0808	
% % % % muv	[1/30,1/14]	0.050985            0.052015                0.0515
% % % % tetathc     6                       5.94                    6.06	
% % % % tetanthc	6                       5.94                    6.06	
% % % % gammaE	[10^3,10^6]                 99000                   101000              100000
% % % % episilonthc     0.7                 0.693                   0.707	
% % % % episilonnthc	0.7                 0.693                   0.707	
% % % % muE     [0.2,0.4]                   0.297                   0.303               0.3
% % % % gammaL          500                 495                     505	
% % % % lamdathc        0.05                0.0495                  0.0505	
% % % % lambdanthc      0.05                0.0495                  0.0505	
% % % % muL     [0.2,0.4]                   0.297                   0.303               0.3
% % % % muP             0.4                 0.396                   0.404	
% % % % C0              0.01                0.0099                  0.0101	
% % % % kF              0.8                 0.792                   0.808	
% % % % r               0.25                0.2475                  0.2525	
% % % % nvnh	[1/3,3]                     1.584                   1.616               1.6
% % % % nv              1                   0.99                	1.01	
% % % % 				
% % % % condiciones iniciales				
% % % % sh(0)	1			
% % % % sv(0)	0.95			
% % % % Iv(0)	0.05			
% % % % 	ceros en las otras variables			

%%%%%3. Evaluar el modelo con los parametros para cada rango se tomaran 50
%%%%%   muestras para cada rango
%%%%%En cada iteracion se calcula un error cuadratico medio entre  la
%%%%%referencia y salida del sistema
%%%%%al finalizar se calula sensibilidad usandos (&M/M)/(&par/par)=RMSD/((medRef-ref)/ref)


function r=resolv5
format long

    %%condiciones del RungeKutta
    tf= 19.5;%%%20 iteraciones iniciando en 0 hasta 19.5 
             %%%(20-0.5) cada iteracion vale 0.5
    h=0.5; % tamaño de paso
    
    %Archivo de salida
    fi = fopen('sensitivity.txt','w+');
    
   
    % % % % condiciones iniciales
        c(1)=0;    %tiempo inicial
        c(2)=1;    %sh(0)
        c(3)=0;    %Ih
        c(4)=0.95; %sv(0)			
        c(5)=0.05; %Iv(0)
        c(6)=0;    %Ethc
        c(7)=0;    %Enthc
        c(8)=0;    %Lthc
        c(9)=0;    %Lnthc
        c(10)=0;   %Pthc
        c(11)=0;   %PnthC
    
    
   

    m=10; %%número de ecuaciones
    
    
    %%arreglos para almacenar los resultados
    Y=[];
    fh=[];
    fe=[];
   
    
    
    %%opcion para buscar sensibilidad en los parametros
    op=24;
   
    
    switch op
        case 1
            refX=0.0000365296803652968; % muh
        case 2 
            refX=1;    % alpha
        case 3
            refX=0.42; % betavh
        case 4
            refX=0.42; % betahv
        case 5
            refX=0.1428;   % ro 
        case 6
            refX=0.001;     % delta
        case 7 
            refX=0.08;      % gammathc
        case 8 
            refX=0.08;      % gammanthc
        case 9
           refX= 0.0515;    % muv 0.0515
        case 10
            refX=6;         % tetathc
        case 11
            refX=6;         % tetanthc
        case 12
            refX=100000;     % gammaE
        case 13
            refX=0.7;      % episilonthc  
        case 14
            refX=0.7;        % episilonnthc
        case 15
            refX= 0.3;      % muE
        case 16
            refX=500;        % gammaL
        case 17
            refX=0.05;       % lamdathc 
        case 18
            refX=0.05;       % lambdanthc    
        case 19
            refX= 0.3;       % muL 
        case 20
            refX= 0.4;      % muP  
        case 21
            refX=0.01;      % C0
        case 22
            refX=0.8;      % kF 
        case 23
            refX=0.25 ;   % r   
        case 24
            refX=1.6;    % nvnh
        case 25
            refX= 1;      % nv 
    end
    
    
            
    %%funcion de referencia, parametros normales  
    resulA=rungeKutta5(tf,c,h,m,op,refX),         
    
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
        c(2)=1;    %sh(0)
        c(3)=0;    %Ih
        c(4)=0.95; %sv(0)			
        c(5)=0.05; %Iv(0)
        c(6)=0;    %Ethc
        c(7)=0;    %Enthc
        c(8)=0;    %Lthc
        c(9)=0;    %Lnthc
        c(10)=0;   %Pthc
        c(11)=0;   %PnthC

    m=10; %%número de ecuaciones
    
    %arreglos para almacenar los resultados
    Y=[];
    resulB=[];
    
    resulB=rungeKutta5(tf,c,h,m,op,X(i)),
     
     
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