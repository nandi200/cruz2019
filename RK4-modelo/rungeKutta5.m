%%encuentra el valor de Y para cada valor de X
%%usando un tamaño de paso h y con valoles
%%iniciales y0 y x0
              
function Y=rungeKutta5(tf,c,h,m,op,X)



        t0=c(1);    %tiempo inicial
        w(1)=c(2);    %sh(0)
        w(2)=c(3);    %Ih
        w(3)=c(4); %sv(0)			
        w(4)=c(5); %Iv(0)
        w(5)=c(6);    %Ethc
        w(6)=c(7);    %Enthc
        w(7)=c(8);    %Lthc
        w(8)=c(9);    %Lnthc
        w(9)=c(10);   %Pthc
        w(10)=c(11);   %PnthC

%%vector respuesta
Y=[];
p=0;
T=[];
    for t=t0:h:tf
       
            for i=1:m
                k1(i)=h*f5(w(1),w(2),w(3),w(4),w(5),w(6),w(7),w(8),w(9),w(10),t,i,op,X);
            end
            
            for i=1:m
                x1=w(1)+k1(1)*0.5;
                x2=w(2)+k1(2)*0.5;
                x3=w(3)+k1(3)*0.5;
                x4=w(4)+k1(4)*0.5;
                x5=w(5)+k1(5)*0.5;
                x6=w(6)+k1(6)*0.5;
                x7=w(7)+k1(7)*0.5;
                x8=w(8)+k1(8)*0.5;
                x9=w(9)+k1(9)*0.5;
                x10=w(10)+k1(10)*0.5;
                k2(i)=h*f5(x1,x2,x3,x4,x5,x6,x7,x8,x9,10,t+h*0.5,i,op,X);
            end
               
            for i=1:m  
                x1=w(1)+k2(1)*0.5;
                x2=w(2)+k2(2)*0.5;
                x3=w(3)+k2(3)*0.5;
                x4=w(4)+k2(4)*0.5;
                x5=w(5)+k2(5)*0.5;
                x6=w(6)+k2(6)*0.5;
                x7=w(7)+k2(7)*0.5;
                x8=w(8)+k2(8)*0.5;
                x9=w(9)+k2(9)*0.5;
                x10=w(10)+k2(10)*0.5;
                
                k3(i)=h*f5(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,t+h*0.5,i,op,X);
            end
                
            for i=1:m    
                x1=w(1)+k3(1);
                x2=w(2)+k3(2);
                x3=w(3)+k3(3);
                x4=w(4)+k3(4);
                x5=w(5)+k3(5);
                x6=w(6)+k3(6);
                x7=w(7)+k3(7);
                x8=w(8)+k3(8);
                x9=w(9)+k3(9);
                x10=w(10)+k3(10);
                k4(i)=h*f5(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,t+h,i,op,X);
            end
               
            
            for i=1:m
                w(i)=w(i)+(k1(i)+2*k2(i)+2*k3(i)+k4(i))/6,
            end
                        
                
            p=p+1;
         Y=[Y;w(1)]; %%sobre la primera ecuación
         T=[T;p];
    
    end
    hold on;
    plot(T,Y(:,1),'b');
   %% plot(T,Y(:,2),'r');
end