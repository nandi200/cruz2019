
function sol=fEDO(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,t,m,op,X)

           %%%parametros del modelo
            muh=0.000037530493; 
            a=0.65;   
            bvh=0.425; 
            bhv=0.425; 
            ro=0.175;   
            muv=0.07916666667;    
            tthc=3.5;  
            tnthc=3.5; 
            emthc=50050;  
            emnthc=50050; 
            ethc=0.7;  
            enthc=0.7;
            muE=0.3;  
            lmthc=250250;    
            lmnthc=250250;
            lthc=0.5;    
            lnthc=0.5;     
            muL=0.3;   
            muP=0.4;  %%no variar
            zi=0.5;   %%no variar
            fi=pi/2;  %%no variar
            gv=0.115;             
            gthc=gv(1+zi*cos(((2*pi)/365)*t+fi));
            gnthc=gv(1+zi*cos(((2*pi)/365)*t+fi));
            c0=0.001;  %%no variar
            kF=0.456;
            r=0.014;
            temp1=(kF*c0*exp(r*t))/(kF+c0*(exp(r*t)-1));
            temp2=((0.571-kF)*c0*exp(r*t))/((0.571-kF)+c0*(exp(r*t)-1));
            Nh=temp1+temp2;

switch op
    case 1
        muh=X;
    case 2
        a=X;
    case 3 
        bvh=X;
    case 4
        bhv=X;
    case 5
        ro=X;
    case 6
        muv=X;
    case 7
        tthc=X;
    case 8
        tnthc=X;
    case 9
       muv=X;
    case 10
       tthc=X; 
    case 11   
       tnthc=X; 
    case 12
       emthc=X;    
    case 13
       emnthc=X;     
    case 14
       ethc=X;      
    case 15
       enthc=X;   
    case 16
       muE=X;  
    case 17
       lmthc=X;    
    case 18
       lmnthc=X;      
    case 19
       lthc=X;  
    case 20
       lnthc=X;    
    case 21
        muL=X;     
    case 22
        gv=X;    
    case 23
        kF=X;  
    case 24
        r=X;      
        
end
aa=exp(r*t);
bb=exp((r*t)-1);
Cf=(kF*C0*aa)/(kF+C0*(bb));
Cm=((1-kF)*C0*exp(r*t))/((1-kF)+C0*(exp(r*t)-1));

Ct=Cf+Cm;

A=0;
    switch m
    
        case 1
            sol=muh*(1-y1)-a*bvh*(1-Ct)*nvnh*y4*y1;
    
        case 2
            sol= a*bvh*(1-Ct)*nvnh*y4*y1-(ro+d+muh)*y2;
        case 3
            sol=(gthc/nv)*y9+(gnthc/nv)*y10-a*bhv*(1-Ct)*y2*y3-(muv-Ct)*y3;
        case 4
            sol=a*bhv*(1-Ct)*y2*y3-(muv+Ct)*y4;
        case 5
            sol=(tthc*(1-((y5+y6)/gE))*nv*(y3+y4))-ethc*y5-(muE+Ct)*y5;
        case 6    
            sol=(tnthc*(1-((y5+y6)/gE))*nv*(y3+y4))-enthc*y6-(muE+Ct)*y6;
        case 7
            sol=ethc*(1-((y7+y4)/gL))*y5-lthc*y7-(muL+A+Ct)*y7;
        case 8
            sol=enthc*(1-((y7+y4)/gL))*y6-lnthc*y8-(muL+Ct)*y8;
        case 9
            sol=lthc*y7-gthc*y9-(muP+Ct)*y9;
        case 10
            sol=lnthc*y8-gnthc*y10-(muP+Ct)*y10;
            
    end

end