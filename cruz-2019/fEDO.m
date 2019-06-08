
function sol=fEDO(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,t,m,op,X)

           %%%parametros del modelo
            muh=; 
            a=;   
            bvh=; 
            bhv=; 
            ro=;   
            d=;   
            gthc=;   
            gnthc=;   
            muv=;    
            tthc=;  
            tnthc=; 
            gE=;   
            ethc=;  
            enthc=; 
            muE=;  
            gL=;    
            lthc=;    
            lnthc=;     
            muL=;   
            muP=;   
            C0=;  
            kF=;  
            r=; 
            nvnh=;  
            nv=;  
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
        d=X;
    case 7
        gthc=X;
    case 8
        gnthc=X;
    case 9
       muv=X;
    case 10
       tthc=X; 
    case 11   
       tnthc=X; 
    case 12
       gE=X;    
    case 13
       ethc=X;     
    case 14
       enthc=X;      
    case 15
            muE=X;   
    case 16
            gL=X;  
    case 17
            lthc=X;    
    case 18
            lnthc=X;      
    case 19
            muL=X;  
    case 20
            muP=X;    
    case 21
            C0=X;     
    case 22
            kF=X;    
    case 23
            r=X;  
    case 24
            nvnh=X;  
    case 25
            nv=X;    
        
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