
function sol=f4(x,y,t,m,op,X)

a=-0.1;
b=0.1;
c=0.2;
switch op
    case 1
        a=X;
    case 2
        b=X;
    case 3 
        c=X;     
end

    switch m
    
        case 1
            sol=a*x;
    
        case 2
            sol= b*x-c*y;
    end

end