%% for binary files read in dec data (singed)
% input: binary file with wordlengh
% output: data contained by variable(matrix) 
function dec_data = file_bin2dec_signed( file_name, bin_length)
    arguments(Input)
        file_name (:,:) char
        bin_length
    end
    pattern_number = floor(length(file_name)/ (bin_length+1));
    for i = 1:pattern_number
        for j = 1:bin_length
            temp(j) = file_name(((bin_length+1)*(i-1))+j);
        end
        dec_temp = bin2dec(temp);
        if(dec_temp >= (2^(bin_length-1)))
            dec_data(i) = dec_temp - (2^(bin_length));
        else 
            dec_data(i) = dec_temp;
        end
    end
end