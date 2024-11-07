%% verilog result 

clear
clc

format long


% file_x = fopen('./xin.txt','w');
% file_x_bin = fopen('./xin_bin.txt','w');
% 
% file_y = fopen('./yin.txt','w');
% file_y_bin = fopen('.yin_bin.txt','w');
% 
% file_ele_angle = fopen('./elementry_ang.txt','w');
% file_ele_angle_bin = fopen('./elementry_ang_bin.txt','w');

file_xout = fopen('./xout.txt','w');
file_xout_bin = fopen('./xout_bin.txt','w');

file_yout = fopen('./yout.txt','w');
file_yout_bin = fopen('.yout_bin.txt','w');

file_ang = fopen('./ang_out.txt','w');
file_ang_bin = fopen('./ang_out_bin.txt','w');

s_n = ones(1,35);       % s1.14 (16 bits)
temp_2 = zeros(1,35);

for i = 1:35
    temp_2(1,i) = 1/(sqrt(1+2^(-2*(i-1))));
end
for i = 1:35
    for j= 1:i
        s_n(1,i) = s_n(1,i)*(temp_2(1,j));
    end
end

trun_sn = 2^12 * truncation(s_n, 12);

sn_11 = s_n(1,11);

bin_sn = dec2bin(trun_sn);

the_sn = bin_sn(11,:);
num_sn = zeros(1, length(the_sn)) ;

bin2dec(the_sn)
CSD_in = append('00',the_sn);
CSD_sn = CSD(CSD_in);



beta = mod(5,4)+1;
alpha = zeros(1,11);

for i = 0:10
    alpha(1,i+1) = (5*i+beta)*pi/25;
end

xin = zeros(1,11);
yin = zeros(1,11);

xout = zeros(1,11);
sxout = zeros(11,1);
yout = zeros(1,11);
out_angle =zeros(1,11);
fix_ele_angle = zeros(20,20);
float_ele_angle = zeros(20,20);


ang = zeros(1,20);
num = zeros(11,1);

for i = 0:10
    xin(1,i+1) = abs(truncation(sin(alpha(1,i+1)), 12));
    xin_orgrin(1,i+1) = truncation(sin(alpha(1,i+1)), 12);
    yin(1,i+1) = truncation(cos(alpha(1,i+1)), 12);
    float_ang(1,i+1) = atan(yin(1,i+1)/xin(1,i+1));
end


% for k = 1:11
%     nbytes_x = fprintf(file_x,'xin[%2d] = %d;\n' , k-1, (2^12)*xin_orgrin(1,k));
%     nbytes_x2 = fprintf(file_x_bin,'xin[%2d] = %s;\n' , k-1, dec2bin((2^12)*xin_orgrin(1,k),14));
% 
%     nbytes_y = fprintf(file_y,'yin[%2d] = %d;\n' , k-1, (2^12)*yin(1,k));
%     nbytes_y2 = fprintf(file_y_bin,'yin[%2d] = %s;\n' , k-1, dec2bin((2^12)*yin(1,k),14));
% end


    for i =1:11
        num(i,1) = 0;
        cnt = 0;
        mu_i = 0;
        x0 = xin(1,i);
        y0 = yin(1,i);
        ang0 = 0;
    
          for m = 1:11
            mu_i = -sign(y0);
    
            x1 = x0;
            y1 = y0;
            ang1 = ang0;
    
            x0 = truncation(x1 - truncation(mu_i*(2^(-cnt))*y1, 12),12);
            y0 = truncation(truncation( mu_i*(2^(-cnt))*x1, 12) + y1,12);
            ang0 = truncation(ang1 - truncation(mu_i*atan(2^(-cnt)),13),13);
            ele_test_angle(i,cnt+1) = 2^13 *truncation(mu_i*atan(2^(-cnt)),12);

            cnt = cnt+1;
            num(i,1) = num(i,1)+1;

    
            xout(i,cnt) = x0;
            yout(i,cnt) = y0;
            ang(i,cnt) = ang0;

          end

    end



for k = 1:11
    nbytes_xout = fprintf(file_xout,'xin[%2d] = %d;\n' , k-1, (2^12)*xout(k,11));
    nbytes_xout2 = fprintf(file_xout_bin,'xin[%2d] = %s;\n' , k-1, dec2bin((2^12)*xin_orgrin(1,k),14));

    nbytes_yout = fprintf(file_yout,'yout[%2d] = %d;\n' , k-1, (2^12)*yout(k,11));
    nbytes_yout2 = fprintf(file_yout_bin,'yout[%2d] = %s;\n' , k-1, dec2bin((2^12)*yout(k,11),14));

    nbytes_ang_out = fprintf(file_ang,'ang[%2d] = %d;\n' , k-1, ((2^13)*truncation(ang(k,11),13)));
    nbytes_ang_out2 = fprintf(file_ang_bin,'ang[%2d] = %s;\n' , k-1, dec2bin((2^13)*ang(k,11),15));
end

    
for i = 1:20 % fraction part
    for j= 1:20 % N step
        fix_ele_angle(i,j) = truncation(atan(2^(-(j-1))),i);
        float_ele_angle(i,j) = atan(2^(-(j-1)));
    end
end


bin_fix_ele_angle = fix_ele_angle(12,:)*(2^13);

% for k =1:11
%     nbytes = fprintf(file_ele_angle,'N[%2d] = %d;\n' , k-1, bin_fix_ele_angle(1,k));
%     nbytes2 = fprintf(file_ele_angle_bin,'N[%2d] = %s;\n' , k-1, dec2bin(bin_fix_ele_angle(1,k),12));
% end


ST = fclose('all');
