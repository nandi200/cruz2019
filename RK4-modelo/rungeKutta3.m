%%encuentra el valor de Y para cada valor de X
%%usando un tamaño de paso h y con valoles
%%iniciales y0 y x0
              
function Y=rungeKutta3(tf,c,h,m)

t0=c(1);
x=c(2);
y=c(3);

w(1)=x;
w(2)=y;
%%vector respuesta
Y=[];
p=0;
T=[];
    for t=t0:h:tf
       
            for i=1:m
                k1(i)=h*f2(w(1),w(2),t,i);
            end
            
            for i=1:m
                k2(i)=h*f2(w(1)+k1(1)*0.5,w(2)+k1(2)*0.5,t+h*0.5,i);
            end
               
            for i=1:m  
                k3(i)=h*f2(w(1)+k2(1)*0.5,w(2)+k2(2)*0.5,t+h*0.5,i);
            end
                
            for i=1:m    
                k4(i)=h*f2(w(1)+k3(1),w(2)+k3(2),t+h,i);
            end
               
            
            for i=1:m
                w(i)=w(i)+(k1(i)+2*k2(i)+2*k3(i)+k4(i))/6,
            end
         
            p=p+1;
         Y=[Y;w];
         T=[T;p];
    
    end
    hold on;
    plot(T,Y(:,1),'b');
    plot(T,Y(:,2),'r');
end