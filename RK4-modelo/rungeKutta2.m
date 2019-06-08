%%encuentra el valor de Y para cada valor de X
%%usando un tamaño de paso h y con valoles
%%iniciales y0 y x0
              
function rungeKutta2(tf,c,h)

t0=c(1);
x=c(2);
y=c(3);

    for t=t0:h:tf
            k1=h*f(x,y,t);
            l1=h*g(x,y,t);
            k2=h*f(x+k1*0.5,y+l1*0.5,t+h*0.5);
            l2=h*g(x+k1*0.5,y+l1*0.5,t+h*0.5);
            k3=h*f(x+k2*0.5,y+l2*0.5,t+h*0.5);
            l3=h*g(x+k2*0.5,y+l2*0.5,t+h*0.5);
            k4=h*f(x+k3,y+l3,t+h);
            l4=h*g(x+k3,y+l3,t+h);
            
            x=x+(k1+2*k2+2*k3+k4)/6,
            y=y+(l1+2*l2+2*l3+l4)/6,
            
        
    
    end
end