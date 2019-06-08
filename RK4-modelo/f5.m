
function sol=f5(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,t,m,op,X)


            muh=0.0000365296803652968; % muh  1
            a=1;    % alphac  2
            bvh=0.42; % betavh 3
            bhv=0.42; % betahv 4 
            ro=0.1428;   % ro  5
            d=0.001;     % delta  6 
            gthc=0.08;      % gammathc  7
            gnthc=0.08;      % gammanthc  8
            muv= 0.0515;    % muv 0.0515  9
            tthc=6;         % tetathc  10
            tnthc=6;         % tetanthc  11
            gE=100000;     % gammaE 12
            ethc=0.7;      % episilonthc  13 
            enthc=0.7;        % episilonnthc  14
            muE= 0.3;      % muE  15
            gL=500;        % gammaL  16
            lthc=0.05;       % lamdathc  17 
            lnthc=0.05;       % lambdanthc 18   
            muL= 0.3;       % muL 19
            muP= 0.4;      % muP  20
            C0=0.01;      % C0 21
            kF=0.8;      % kF 22
            r=0.25 ;   % r   23
            nvnh=1.6;    % nvnh  24
            nv= 1;      % nv 25
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
       tnthc=X;   % tetanthc  11
    case 12
       gE=X;     % gammaE 12
    case 13
       ethc=X;      % episilonthc  13
    case 14
       enthc=X;        % episilonnthc  14
    case 15
            muE=X;      % muE  15
    case 16
            gL=X;        % gammaL  16
    case 17
            lthc=X;       % lamdathc  17
    case 18
            lnthc=X;       % lambdanthc 18 
    case 19
            muL=X;       % muL 19
    case 20
            muP=X;      % muP  20
    case 21
            C0=X;      % C0 21
    case 22
            kF=X;      % kF 22
    case 23
            r=X;   % r   23
    case 24
            nvnh=X;    % nvnh  24
    case 25
            nv=X;      % nv 25
        
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