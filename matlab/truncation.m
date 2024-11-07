function T_z = truncation ( z , a ) %z = value , a = fraction number
    format long
    if (a == 0) 
        T_z = z;
    else 
        T_z = floor(z*(2^a))/(2^a);
    end
end
