
function sol=f2(x,y,t,m)

    switch m
    
        case 1
            sol=-0.1*x;
    
        case 2
            sol= 0.1*x-0.2*y;
    end

end