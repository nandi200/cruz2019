
function sol=fEDO(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,t,m,op,X)

           %%%parametros del modelo
            muh=0.000037530493; 
            a=0.65;   
            bvh=0.425; 
            bhv=0.425; 
            ro=0.025;   
            muv=0.0775;    
            tthc=3.5;  
            tnthc=3.5; 
            emthc=500500;  
            emnthc=500500; 
            ethc=0.7;  
            enthc=0.7;
            muE=0.3;  
            lmthc=250250;    
            lmnthc=250250;
            lthc=0.5;    
            lnthc=0.5;     
            muL=0.3;   
            muP=0.4;
            A=0.1;
            kmax=0.573;  %%kF
            pf=0.13;
            pm=0.04;
            c0=0.001;  %%no variar
            r=0.015;
            zi=0.5; 
            fi=pi/2;
            g0=0.115;     %%gv      
            gthc=gv(1+zi*cos(((2*pi)/365)*t+fi));
            gnthc=gv(1+zi*cos(((2*pi)/365)*t+fi));
           
           %% kF=0.456;
           
            temp1=(kF*c0*exp(r*t))/(kF+c0*(exp(r*t)-1));
            temp2=((0.571-kF)*c0*exp(r*t))/((0.571-kF)+c0*(exp(r*t)-1));
            Ct=temp1+temp2;
            Nh=834634;

switch op
    case 1
        muh=X;  %%no variar
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
       emthc=X;
    case 10
       emnthc=X; 
    case 11   
       ethc=X; 
    case 12
       enthc=X;    
    case 13
       muE=X;     
    case 14
       lmthc=X;      
    case 15
       lmnthc=X;   
    case 16
       lthc=X;  
    case 17
       lnthc=X;    
    case 18
       muL=X;      
    case 19   %%no variar
       muP=X;  
    case 20   %%no variar
       zi=X;    
    case 21  %%no variar
       fi=X;     
    case 22
       gv=X;    
    %%falta
    %%gthc
    %%gnthc
    case 23
       c0=X;  
    case 24
       kF=X;       
    case 25
       r=X;
    case 26
       Nh=x;
        
end
%aa=exp(r*t);
%bb=exp((r*t)-1);
%Cf=(kF*C0*aa)/(kF+C0*(bb));
%Cm=((0.571-kF)*C0*exp(r*t))/((0.571-kF)+C0*(exp(r*t)-1));
%Ct=Cf+Cm;

%A=0;
    switch m
    
        case 1  %%Sh
            sol=muh*(Nh-y1)-a*(bvh/Nh)*(1-Ct)*y4*y1;
    
        case 2  %%Ih
            sol= a*(bvh/Nh)*(1-Ct)*y4*y1-(ro+muh)*y2;
            
        case 3 %%Sv
            sol=(gthc*y9)+(gnthc*y10)-a*(bhv/Nh)*(1-Ct)*y2*y3-(muv+(muv*Ct))*y3;
            
               %%Iv
        case 4
            sol=a*(bhv/Nh)*(1-Ct)*y2*y3-(muv+(muv*Ct))*y4;
            
               %%Ethc
        case 5
            sol=(tthc*(1-(y5/emthc))*(y3+y4))-ethc*y5-(muE+(muE*Ct))*y5;
            
              %%Enthc
        case 6    
            sol=(tnthc*(1-(y6/emnthc))*(y3+y4))-enthc*y6-(muE+(muE*Ct))*y6;
            
                %%Lthc
        case 7
            sol=ethc*(1-(y7/lmthc))*y5-lthc*y7-(muL+A+(muL*Ct))*y7;
            
        case 8  %%Lnthc
            sol=enthc*(1-(y7/lmthc))*y6-lnthc*y8-(muL+(muL*Ct))*y8;
                        
        case 9  %%Pthc
            sol=lthc*y7-gthc*y9-(muP+(muP*Ct))*y9;
            
        case 10  %%Pnthc
            sol=lnthc*y8-gnthc*y10-(muP+(muP*Ct))*y10;
            
        case 11 %%hgthc=hgnthc
            kk=(sin(fi)*cos((2*pi*t)/365))+(sin(2*pi*t/365)*cos(fi));
            sol=muv+(365/*pi*t)*zi*muv*(kk-sin(fi));
            
        case 12 %%hC
            sol=1/t(((kF/r)*log(1+(c0*(exp(r*t)-1)/kF)))+((0.571-kF)/r)*log(1+((c0*(exp(r*t)-1))/(0.571-kF))); )
            
    end
    
            aa=(enthc*emnthc+(lnthc+muL*y12+muL)*lmnthc)*(ethc*emthc+(lthc+muL*y12+A+muL)*lmthc);
            bb=lmnthc*emnthc*lmthc*emthc;
            A1=aa/bb*(muv*y12+muv)*tthc*tnthc;   %%deben ser mayor a 0;
    
            cc=(lnthc+muL*y12+muL)*(enthc+muE*y12+muE)*(muv*y12+muv);
            dd=(y11*lmthc*enthc*tnthc)/(y12+muv*y11+muv);
            ee=((ethc/lmthc)+(lthc+muL*y12+A+muL)/emthc)*tthc;
            ff=(lthc+muL*y12+A+muL)*(ethc+muE*y12+muE)*(muv*y12+muv);
            gg=(y11*lthc*ethc*tthc)/(y11+muP*y12+muP);
            hh=((enthc/lmnthc)+((lnthc+muL*y12+muL)/emnthc))*tnthc;  %%debe ser mayor a 0;
            
            A2=(cc-dd)*(ee)+(ff-gg)*(hh);
            
            ii=(muv*y12+muv)*(lthc+muL*y12+A+muL)*(ethc+muE*y12+muE)*(y11+muP*y12*muP);
            RthcM=y11*lthc*ethc*tthc/(ii);
            jj=(muv*y12+muv)*(lnthc+muL*y12+muL)*(enthc+muE*y12+muE)*(y11+muP*y12+muP);
            RnthcM=y11*lnthc*enthc*tnthc;
            Rm=RthcM+RnthcM/jj;
            A3=(muv*y12+muv)*(lthc+muL*y12+A+muL)*(ethc+muE*y12+muE)*(lnthc+muL*y12+muL)*(enthc+muE*y12+muE)*(1-Rm);
               
            
           M=-A2+pown
       

end