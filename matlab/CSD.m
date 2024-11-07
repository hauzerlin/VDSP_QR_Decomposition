function c_ans = CSD (input_b) 
    format long

    length_b = length(input_b);

    s = zeros(1,length_b);
    g = zeros(1,length_b+1);
    one_min_two_b = zeros(1,length_b);
    num_sn = zeros(1,length_b+2);
    

    for i=1:length_b
        if(input_b(length_b+1-i)=='1')
            num_sn(i) = 1;
        else
            num_sn(i) = 0;
        end
    end
    num_sn(length_b+1) = 0;
    num_sn = circshift(num_sn,1);

    for i = 2:length_b+1
        s(i-1) = xor(num_sn(i),num_sn(i-1));
    end

    for i = 2:length_b+1
        g(i) = and(~g(i-1),s(i-1));
    end

       for i =3: length_b+2
        if(num_sn(i)==1)
            one_min_two_b(i-2) = -1;
        else
            one_min_two_b(i-2) = 1;
        end
       end

    for i=1:length_b
        c(i) = one_min_two_b(i)*g(i+1);
    end

    c_ans = flip(c);
end
