function [B, phi] = extQfromDCM(BN)

B(1) = 1/2*sqrt(1+trace(BN));
B(2) = 1/2*(1+2*BN(1,1)-trace(BN));
B(3) = 1/2*sqrt(1+2*BN(2,2)-trace(BN));
B(4) = 1/2*sqrt(1+2*BN(3,3)-trace(BN));
[Bmax, I] = max(B);
    if I == 1
        B(2) = (BN(2,3)-BN(3,2))/(4*B(1));
        B(3) = (BN(3,1)-BN(1,3))/(4*B(1));
        B(4) = (BN(1,2)-BN(2,1))/(4*B(1));
    elseif I==2
        B(1) = (BN(2,3)-BN(3,2))/(4*B(2));
        B(3) = (BN(1,2)+BN(2,1))/(4*B(2));
        B(4) = (BN(3,1)+BN(1,3))/(4*B(2));
    elseif I==3
        B(1) = (BN(3,1)-BN(1,3))/(4*B(3));
        B(2) = (BN(1,2)+BN(2,1))/(4*B(3));
        B(4) = (BN(2,3)+BN(3,2))/(4*B(3));
    else
        B(1) = (BN(1,2)-BN(2,1))/(4*B(4));
        B(2) = (BN(3,1)+BN(1,3))/(4*B(4));
        B(3) = (BN(2,3)+BN(3,2))/(4*B(4));
    end

    phi  = 2*acos(B(1));
    if phi > pi
        phi = phi - 2*pi;
    end
    
end