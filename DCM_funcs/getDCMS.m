function [BN, dcms] = getDCMS(sequence,theta)
    dcms = struct('DCMS',{});
    for i = 1:3
       if sequence(i) == 1
           dcms(i).DCMS = M1(theta(i));
       elseif sequence(i) == 2
           dcms(i).DCMS = M2(theta(i));   
       else
           dcms(i).DCMS = M3(theta(i));
       end
    end
    
    BN = dcms(3).DCMS*dcms(2).DCMS*dcms(1).DCMS;
end




