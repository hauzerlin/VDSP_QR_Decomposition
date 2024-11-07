%% For normal CORDIC function
clear;clc;

% stage of Rotation
Micro_Rotation = 11;
Angle_Accumulation = 11;
stage = 11;

% Wordlength
X_WL = 12; % S1.12;
Y_WL = 12; % S1.12;
element_angle_WL = 13; % s1.13;
Sn_WL = 12; % S1.12;
CSD_WL = 12;

% Parameter
data_pattern = 11;


beta = mod(5,4)+1; % student ID = 111521035, I = the last number of student ID = 5
i = 0:data_pattern-1;
alpha = (5*i+beta)*pi/25;


N = 35;
s_n = ones(1,N);
temp_2 = zeros(1,N);

% S(N) = 1/Π(1+2^(-2i))
for i = 1:N
    temp_2(1,i) = 1/(sqrt(1+2^(-2*(i-1)))); 
end
for i = 1:N
    for j= 1:i
        s_n(1,i) = s_n(1,i)*(temp_2(1,j));
    end
end
s_n = truncation(s_n,Sn_WL);

% input data
xin = zeros(1,data_pattern); 
yin = zeros(1,data_pattern);
for i = 0:data_pattern-1
    xin(1,i+1) = abs(truncation(sin(alpha(1,i+1)), X_WL));
    yin(1,i+1) = truncation(cos(alpha(1,i+1)), Y_WL);
end

% elementary angle
for idk = 0:Micro_Rotation-1
    elementary_angle(idk+1) = truncation(atan(1/(2^idk)),element_angle_WL);
end

xout = zeros(data_pattern,data_pattern); % save x states
yout = zeros(data_pattern,data_pattern); % save y states
mag_error = zeros(data_pattern,stage); % save magnitude errors

ang = zeros(1,stage); % save angle states
num = zeros(data_pattern,1); % to count the times that roate the angle

for i =1:data_pattern
    
    num(i,1) = 0;
    cnt = 0;
    mu_i = 0;
    x0 = xin(1,i);
    y0 = yin(1,i);
    ang0 = 0;

    for k= 1:stage
        % to rotate the angle
        mu_i = -sign(y0);

        x1 = x0;
        y1 = y0;
        ang1 = ang0;

        % X and Y are quantized into S1.12
        x0 = x1 - truncation(mu_i*(2^(-cnt))*y1, 12);
        y0 = truncation( mu_i*(2^(-cnt))*x1, 12) + y1;
        ang0 = ang1 - mu_i*truncation(atan(2^(-cnt)),element_angle_WL);

        cnt = cnt+1;
        num(i,1) = num(i,1)+1;
        
        % to save the rotation result
        xout(i,cnt) = x0;
        yout(i,cnt) = y0;
        ang(i,cnt) = ang0;
    end
end

for idx = 1: size(xout,2)
x_mapping(idx,:) = xout(idx,:) .* s_n(size(xout,1));
end



%% test CSD for autocontrolling

% clc;clear;

N = 35;
s_n = ones(1,N);
temp_2 = zeros(1,N);

% S(N) = 1/Π(1+2^(-2i))
for i = 1:N
    temp_2(1,i) = 1/(sqrt(1+2^(-2*(i-1)))); 
end
for i = 1:N
    for j= 1:i
        s_n(1,i) = s_n(1,i)*(temp_2(1,j));
    end
end

one_of_s_n = 1./s_n;
% 0.607252935008881
CSD_data =  CSD(['0',dec2bin((2^13)*truncation(0.607252935008881,13))]);
SN_length = 13;
x_CSD = 0;
for idx =2:length(CSD_data)
                if(CSD_data(idx)==1)
                   x_CSD = x_CSD + ...
                       truncation(1* ...
                       truncation(2^-(idx-1),SN_length),SN_length); 
               elseif(CSD_data(idx)==-1)
                   x_CSD = x_CSD - ...
                       truncation(1* ...
                       truncation(2^-(idx-1),SN_length),SN_length); 
                end
end

%% test vectoring mode with pattern grenrated by wayne
clc;clear;

% xin_file_wayne = fopen('txt/input_pat_x.txt');
% xin_wayne = fscanf(xin_file_wayne, '%15c');

pattern = 10000;
input_length =15;

xin_wayne_bin = fileread('txt/input_pat_X.txt');
yin_wayne_bin = fileread('txt/input_pat_Y.txt');

xin_wayne = file_bin2dec_signed(xin_wayne_bin, 14);
yin_wayne = file_bin2dec_signed(yin_wayne_bin, 14);

xout_wayne_bin = fileread('txt/output_pat_X.txt');
yout_wayne_bin = fileread('txt/output_pat_Y.txt');

xout_wayne = file_bin2dec_signed(xout_wayne_bin, 14);
yout_wayne = file_bin2dec_signed(yout_wayne_bin, 14);

fclose('all');
xin = xin_wayne ./ (2^12);
yin = yin_wayne ./ (2^12);
n = pattern; 
%

for i =1:n % run 11 inputs
  [rad_data(i), x_data(i), y_data(i), mapping_flag(i)] = ...
      CORDIC( ...
      xin(i), yin(i), 0, stage = 11, ...
      input_length = 12, input_angle_length = 13, ...
      element_angle_length=13 ,SN_length=13 ...
      );
end
angle_data = rad2deg(rad_data);
for i =1:n % run 11 inputs
    golden_angle(i) = atan(yin(i)/xin(i));
end
% turning back
for i = 1:n
    [ angle_test(i), x_test(i), y_test(i) ] = ...
        CORDIC(x_data(i), y_data(i), rad_data(i), mapping_flag(i), 1, ...
      stage= 11, input_length = 12, input_angle_length = 13, ...
      element_angle_length=13 ,SN_length=13 ...
      );
        % CORDIC_180(x_data(i), y_data(i), (10^-100), 0, 1);
end

xout_dec = (2^12) * x_data;
yout_dec = (2^12) * y_data;
angle_dec = (2^13) * rad_data;

ang_error = sqrt(mean((abs(rad_data - golden_angle)).^2));
x_mag_error = sqrt(mean((abs(xin - x_test)).^2));
y_mag_error = sqrt(mean((abs(yin - y_test)).^2));


%% test vectoring mode with pattern grenrated by wayne(offical spec)
clc;clear;

% xin_file_wayne = fopen('txt/input_pat_x.txt');
% xin_wayne = fscanf(xin_file_wayne, '%15c');

pattern = 10000;
input_length =15;

xin_wayne_bin = fileread('txt/input_pat_X.txt');
yin_wayne_bin = fileread('txt/input_pat_Y.txt');

xin_wayne = file_bin2dec_signed(xin_wayne_bin, 14);
yin_wayne = file_bin2dec_signed(yin_wayne_bin, 14);

xout_wayne_bin = fileread('txt/output_pat_X.txt');
yout_wayne_bin = fileread('txt/output_pat_Y.txt');

xout_wayne = file_bin2dec_signed(xout_wayne_bin, 14);
yout_wayne = file_bin2dec_signed(yout_wayne_bin, 14);

fclose('all');
xin = xin_wayne ./ (2^14);
yin = yin_wayne ./ (2^14);
n = pattern; 
%

for i =1:n % run 11 inputs
  [rad_data(i), x_data(i), y_data(i), mapping_flag(i)] = ...
      CORDIC( ...
      xin(i), yin(i), 0, stage = 13, ...
      input_length = 14, input_angle_length = 14, ...
      element_angle_length=14 ,SN_length=14 ...
      );
end
angle_data = rad2deg(rad_data);
for i =1:n % run 11 inputs
    golden_angle(i) = atan(yin(i)/xin(i));
end
% turning back
for i = 1:n
    [ angle_test(i), x_test(i), y_test(i) ] = ...
        CORDIC(x_data(i), y_data(i), rad_data(i), mapping_flag(i), 1, ...
      stage= 13, input_length = 14, input_angle_length = 14, ...
      element_angle_length=14 ,SN_length=14 ...
      );
        % CORDIC_180(x_data(i), y_data(i), (10^-100), 0, 1);
end

xout_dec = (2^14) * x_data;
yout_dec = (2^14) * y_data;
angle_dec = (2^14) * rad_data;

ang_error = sqrt(mean((abs(rad_data - golden_angle)).^2));
x_mag_error = sqrt(mean((abs(xin - x_test)).^2));
y_mag_error = sqrt(mean((abs(yin - y_test)).^2));
%% test rotation mode with angle 0
clc;clear;
[ angle_test, x_test, y_test, flag_test ] = ...
CORDIC(-11926/(2^14), 11922/(2^14), 0/(2^14) , 0, 0, ...
stage= 13, input_length = 14, input_angle_length = 14, ...
element_angle_length=14 ,SN_length=14 ...
);
x_test * (2^14)
y_test * (2^14)
angle_test * (2^14)
flag_test
%% elementary angle 
% SN_wordlenght = S1.14
% stage = 13
clc;clear

N = 20;
SN_wordlength = 14;
stage = 13;

elementary_angles = zeros(1,N+1);
for i = 1:N+1
    elementary_angles(1,i) =  ...
        truncation(atan(1/(2^(i-1))),SN_wordlength);
end
s_n = ones(1,N);
for i = 1:N
    temp_2(1,i) = 1/(sqrt(1+2^(-2*(i-1)))); 
end
for i = 1:N
    for j= 1:i
        s_n(1,i) = s_n(1,i)*(temp_2(1,j));
    end
end
CSD_data = CSD( ['0', dec2bin((2^SN_wordlength)*truncation(s_n(stage),SN_wordlength))] )

%% output for binary data
file_elementary_angle = fopen("hower\elementary_angle_binary.txt","w");
for idx = 1:stage
    nbyte = fprintf(file_elementary_angle,'%2d\''b%s\n',SN_wordlength,...
        dec2bin((2^SN_wordlength)*elementary_angles(idx),SN_wordlength));
end
fclose('all');
%% function for read binray data into varible (vector)

% function dec_data = file_bin2dec_signed( file_name, bin_length)
%     arguments(Input)
%         file_name (1,:) char
%         bin_length
%     end
%     pattern_number = floor(length(file_name)/ (bin_length+1));
%     for i = 1:pattern_number
%         for j = 1:bin_length
%             temp(j) = file_name(((bin_length+1)*(i-1))+j);
%         end
%         dec_temp = bin2dec(temp);
%         if(dec_temp >= (2^(bin_length-1)))
%             dec_data(i) = dec_temp - (2^(bin_length));
%         else 
%             dec_data(i) = dec_temp;
%         end
%     end
% end

%% function for CORDIC.
% function [angle_out, x_out, y_out, flag_out] = ...
%          CORDIC(x, y, angle_in, flag_in, mode, options)
% % mode 0:vectoring mode(clockwise), mode 1:rotation mode(anticlockwise)
% % offical use for final project
% 
%     arguments(Input)
%         x (1,:) double
%         y (1,:) double
%         angle_in double = 0
%         flag_in {mustBeInteger} = 0
%         mode (1,1) = 0 
%         options.stage (1,1) {mustBePositive} = 200
%         options.input_length (1,1) {mustBeInteger} = 0 %input wordlength
%         options.input_angle_length (1,1) {mustBeInteger} = 0
%         options.element_angle_length (1,1) {mustBeInteger} = 0
%         options.SN_length (1,1) {mustBeInteger} = 0
%     end
%     arguments(Output)
%         angle_out
%         x_out
%         y_out
%         flag_out
%     end
%     format long
% 
%     % S(N) = 1/Π(1+2^(-2i))
%     s_n = ones(1,options.stage+1);
%     temp = zeros(1,options.stage+1);
%     for i = 1:options.stage+1
%         temp(1,i) = 1/(sqrt(1+2^(-2*(i-1)))); 
%     end
%     for i = 1:options.stage+1
%         for j= 1:i
%             s_n(1,i) = s_n(1,i)*(temp(1,j));
%         end
%     end
%     s_n = truncation(s_n,options.SN_length);
% 
%     % elementary angle 
%     elementary_angles = zeros(1,options.stage+1);
%     for i = 1:options.stage+1
%         elementary_angles(1,i) =  ...
%             truncation(atan(1/(2^(i-1))), options.element_angle_length);
%     end
% 
%     % initial options.stage
%     x_data = zeros(1,options.stage+1);
%     y_data = zeros(1,options.stage+1);
%     angle_data = zeros(1,options.stage+1);
%     init_angle = truncation(angle_in,options.input_angle_length);
% 
%         if(mode==0)
%             if(x<0)
%                 mapping_flag = 1;
%             else
%                 mapping_flag = 0;
%             end
%             x_data(1) = abs(truncation(x,options.input_length));
%             y_data(1) = truncation(y,options.input_length);
%             init_angle = 0;
%             angle_data(1) = init_angle;
% 
%         elseif (mode==1)
%             if(flag_in==1)
%                 x_data(1) = -truncation(x,options.input_length);
%                 y_data(1) = -truncation(y,options.input_length);
%             else
%                 x_data(1) = truncation(x,options.input_length);
%                 y_data(1) = truncation(y,options.input_length);
%             end
%             angle_data(1) = init_angle;
%             mapping_flag = flag_in;
%         else 
%             error("mode is wrong");
%         end
% 
%     % caculate stage
%     for idx = 1:options.stage
%         if (mode==0) 
%             if(y_data(idx)<0) mu = 1;
%             else mu = -1;
%             end
%             % mu = -sign(y_data(idx));
%         else 
%             if(angle_data(idx)<0) mu=-1;
%             else mu = 1;
%             end
%             % mu = sign(angle_data(idx));
%         end
%         x_data(idx+1) = truncation (x_data(idx) - ...
%                         mu*...
%                         truncation( ...
%                         ( truncation( (2^(-(idx-1)) ),options.input_length)* ...
%                         y_data(idx) ),options.input_length) ,options.input_length);
%         y_data(idx+1) = truncation(mu* ...
%                         truncation( (x_data(idx)* ...
%                         truncation( (2^(-(idx-1)) ),options.input_length) ), ...
%                         options.input_length ) ...
%                         + y_data(idx) ,options.input_length);
%         angle_data(idx+1) = angle_data(idx) - mu*elementary_angles(idx);
%     end
% 
%     % final stage(output): CSD and mapping
%     if(mode==0) % vectroing mode output
%         if(mapping_flag==1)
%             angle_out = -angle_data(options.stage+1);
%         else
%             angle_out = angle_data(options.stage+1);
%         end
%         if(options.SN_length==0)
%             x_out = x_data(options.stage+1) * s_n(options.stage);
%             y_out = y_data(options.stage+1) * s_n(options.stage);
%         else
%             CSD_data = CSD( ['0', dec2bin((2^options.SN_length)* ...
%                 truncation(s_n(options.stage),options.SN_length))] );
%             x_CSD = 0; y_CSD = 0;
%             for idx =2:length(CSD_data)
%                 if(CSD_data(idx)==1)
%                    x_CSD = x_CSD + ...
%                        truncation(x_data(options.stage+1)* ...
%                        truncation(2^-(idx-1),options.input_length),options.input_length); 
%                    y_CSD = y_CSD + ...
%                        truncation(y_data(options.stage+1)* ...
%                        truncation(2^-(idx-1),options.input_length),options.input_length); 
%                 elseif(CSD_data(idx)==-1)
%                    x_CSD = x_CSD - ...
%                        truncation(x_data(options.stage+1)* ...
%                        truncation(2^-(idx-1),options.input_length),options.input_length); 
%                    y_CSD = y_CSD - ...
%                        truncation(y_data(options.stage+1)* ...
%                        truncation(2^-(idx-1),options.input_length),options.input_length); 
%                 end
%             end
%             x_out = truncation(x_CSD,options.input_length);
%             y_out = truncation(y_CSD,options.input_length);
%         end
%         flag_out  = mapping_flag;
% 
%     elseif (mode==1) % rotation mode output
%         angle_out = angle_data(options.stage+1);
% 
%         % x_out = x_data(options.stage+1) * s_n(options.stage);
%         % y_out = y_data(options.stage+1) * s_n(options.stage);
%         if(options.SN_length==0)
%             x_out = x_data(options.stage+1) * s_n(options.stage);
%             y_out = y_data(options.stage+1) * s_n(options.stage);
%         else
%             CSD_data = CSD( ['0', dec2bin((2^options.SN_length)* ...
%                 truncation(s_n(options.stage),options.SN_length))] );
%             x_CSD = 0; y_CSD = 0;
%             for idx =2:length(CSD_data)
%                 if(CSD_data(idx)==1)
%                    x_CSD = x_CSD + ...
%                        truncation(x_data(options.stage+1)* ...
%                        truncation(2^-(idx-1),options.input_length),options.input_length); 
%                    y_CSD = y_CSD + ...
%                        truncation(y_data(options.stage+1)* ...
%                        truncation(2^-(idx-1),options.input_length),options.input_length); 
%                 elseif(CSD_data(idx)==-1)
%                    x_CSD = x_CSD - ...
%                        truncation(x_data(options.stage+1)* ...
%                        truncation(2^-(idx-1),options.input_length),options.input_length); 
%                    y_CSD = y_CSD - ...
%                        truncation(y_data(options.stage+1)* ...
%                        truncation(2^-(idx-1),options.input_length),options.input_length); 
%                 end
%             end
%             x_out = truncation(x_CSD,options.input_length);
%             y_out = truncation(y_CSD,options.input_length);
%         end
% 
%         flag_out  = flag_in;
%     end
% % s_n(options.stage)
% end


%% function for CSD
% function c_ans = CSD (input_b) 
%     format long
% 
%     length_b = length(input_b);
% 
%     s = zeros(1,length_b);
%     g = zeros(1,length_b+1);
%     one_min_two_b = zeros(1,length_b);
%     num_sn = zeros(1,length_b+2);
% 
% 
%     for i=1:length_b
%         if(input_b(length_b+1-i)=='1')
%             num_sn(i) = 1;
%         else
%             num_sn(i) = 0;
%         end
%     end
%     num_sn(length_b+1) = 0;
%     num_sn = circshift(num_sn,1);
% 
%     for i = 2:length_b+1
%         s(i-1) = xor(num_sn(i),num_sn(i-1));
%     end
% 
%     for i = 2:length_b+1
%         g(i) = and(~g(i-1),s(i-1));
%     end
% 
%        for i =3: length_b+2
%         if(num_sn(i)==1)
%             one_min_two_b(i-2) = -1;
%         else
%             one_min_two_b(i-2) = 1;
%         end
%        end
% 
%     for i=1:length_b
%         c(i) = one_min_two_b(i)*g(i+1);
%     end
% 
%     c_ans = flip(c);
% end

%% function for truncation
% 
% function T_z = truncation ( z , a ) %z = value , a = fraction number
%     format long
%     if (a == 0) 
%         T_z = z;
%     else 
%         T_z = floor(z*(2^a))/(2^a);
%     end
% end
