%% procedure1 
% Although we do not perform vector scaling when using CORDIC to obtain the
% phase of a complex number, we still need to know the increase in 
% magnitude after infinite micro-rotations for reserving sufficient dynamic
% range during implementation. Find out the scaling factor 𝑆(𝑁) for ≥20.

clear;clc;
format long

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

% plot(s_n(1:5));
plot(s_n);
xlabel('N'),ylabel('S(N) value');
title('scaling factor S(N)');

%% procedure_2 practice
% Assume that 𝑋=sin (𝛼) , 𝑌=cos (𝛼) , where 𝛼=(5𝑛+𝛽)𝜋/25 for 𝑛=0,1,…,10 
% and 𝛽=𝑚𝑜𝑑(𝐼,4)+1, where 𝐼 is the last digit in your student ID.
% Both 𝑋 and 𝑌 are quantized into 14 bits including the sign bit and the
% 12-bit fractional part. 
% According to Q1, determine the word-length of 𝑋(𝑖) and 𝑌(𝑖) at all the 
% stages if they use the same format 
% (Hint: Consider the possible growth of the input signal.)

clear;clc;
n = 11;

% student ID = 111521035, I = the last number of student ID = 5
beta = mod(5,4)+1;
alpha = zeros(1,n);

i = 0:n-1;
alpha = (5*i+beta)*pi/25;
% for i = 0:10
%     alpha(1,i+1) = (5*i+beta)*pi/25;
% end

xin = zeros(1,n);
yin = zeros(1,n);
xout = zeros(1,n);
yout = zeros(1,n);

test_ang = zeros(1,n);

Stage = 20; % 試試看使用多少wordlength足夠長

ang = zeros(1,Stage);
num = zeros(n,1);
shift_ang = zeros(n,Stage);

for i = 0:n-1
    xin(1,i+1) = abs(truncation(sin(alpha(1,i+1)), 12));
%     xin(1,i+1) = (truncation(sin(alpha(1,i+1)), 12));
    yin(1,i+1) = truncation(cos(alpha(1,i+1)), 12);
%     test_ang(1,i+1) = atan(yin(1,i+1)/xin(1,i+1));
%     rev_angle(1,i+1) = atan (yin(1,i+1)/abs(xin(1,i+1)));
end

for i =1:n % run 11 inputs
    
    % initial 
    num(i,1) = 0;
    cnt = 0;
    mu_i = 0; 
    x0 = xin(1,i);
    y0 = yin(1,i);
    ang0 = 0;

    for j =1 :Stage % run 20 stages
        
        mu_i = -sign(y0); % check if y0>0

        x1 = x0;
        y1 = y0;
        ang1 = ang0;

        x0 = x1 - truncation(mu_i*(2^(-cnt))*y1, 12);
        y0 = truncation( mu_i*(2^(-cnt))*x1, 12) + y1;
        ang0 = ang1 - mu_i*atan(2^(-cnt));

        shift_ang(i,j) = abs(mu_i*atan(2^(-cnt))); % |位移角度|
        cnt = cnt+1;
        num(i,1) = num(i,1)+1;

        xout(i,j) = x0; % x0 每一次旋轉的變化紀錄
        yout(i,j) = y0; % y0 每一次旋轉的變化紀錄
        ang(i,cnt) = ang0; % angle 每一次旋轉的變化紀錄

    end
end



%% procedure_2 procedure

clear;clc;

n = 11; % number of input data is 11

% student ID = 111521035, I = the last number of student ID = 5
beta = mod(5,4)+1;
alpha = zeros(1,11);

i = 0:n-1;
alpha = (5*i+beta)*pi/25;

xin = zeros(1,11);
yin = zeros(1,11);
xout = zeros(1,11);
yout = zeros(1,11);

ang = zeros(1,20);
float_ang = zeros(1,11);
floting_error = zeros(11,20);
num = zeros(11,1);
shift_ang = zeros(11,20);

for i = 0:10
    xin(1,i+1) = abs(truncation(sin(alpha(1,i+1)), 12));
    xin_orgrin(1,i+1) = truncation(sin(alpha(1,i+1)), 12);
    yin(1,i+1) = truncation(cos(alpha(1,i+1)), 12);
    float_ang(1,i+1) = atan(yin(1,i+1)/xin(1,i+1));
end

for k = 1:16 % WL of truncation
    for i =1:11 % 11 inputs

        % initial
        num(i,1) = 0;
        cnt = 0;
        mu_i = 0;
        x0 = xin(1,i);
        y0 = yin(1,i);
        ang0 = 0;

        for j =1 :20  % 20-stage roations
            mu_i = -sign(y0);
    
            x1 = x0;
            y1 = y0;
            ang1 = ang0;

            % X and Y are quantized into S1.12
            x0 = x1 - truncation(mu_i*(2^(-cnt))*y1, 12);
            y0 = truncation((truncation(mu_i*(2^(-cnt))*x1, 12) + y1),k);
            ang0 = ang1 - mu_i*atan(2^(-cnt));

            cnt = cnt+1;
            num(i,1) = num(i,1)+1;

            ang(i,cnt) = ang0;
        end
        
        % RMSE 中的 Square and Error
        floting_error(i,k) = (float_ang(1,i) - ang(i,20))^2;

    end
end

% RMSE 中的 Root and mean
avg_phase_err = sqrt(mean(floting_error));

for avg_cnt = 1:20
    if((avg_phase_err(1,avg_cnt)) <(0.4*2^(-9)))
        break
    end
end


plot(avg_phase_err)
set(gca, 'YScale', 'log')
% title('different word-length of X(i) versus the root mean squared error');
title('different word-length of Y(i) versus the root mean squared error');
xlabel('word-length(bits)'),ylabel('RMSE of output');
% title('The average phase error');
% xlabel('numbers of micro-rotations 𝑁'), ylabel('the average phase errors')
yline(0.4*2^(-9),'-r','0.4*2 ^-^9')

% file_x = fopen('./xin.txt','w');
% file_x_bin = fopen('./xin_bin.txt','w');
% file_y = fopen('./yin.txt','w');

% file_y_bin = fopen('.yin_bin.txt','w');
% for k = 1:11
%     nbytes_x = fprintf(file_x,'xin[%2d] = %d;\n' , k-1, (2^12)*xin_orgrin(1,k));
%     nbytes_x2 = fprintf(file_x_bin,'xin[%2d] = %s;\n' , k-1, dec2bin((2^12)*xin_orgrin(1,k),14));
% 
%     nbytes_y = fprintf(file_y,'yin[%2d] = %d;\n' , k-1, (2^12)*yin(1,k));
%     nbytes_y2 = fprintf(file_y_bin,'yin[%2d] = %s;\n' , k-1, dec2bin((2^12)*yin(1,k),14));
% end

% ST = fclose('all');

%% procedure_3 practice
% 決定 numbers of micro-rotations的級數使得elementary angle的 RMSE 在容忍值內
% average phase error

clear;clc;

n = 11; % number of input data is 11
stage = 20; % max number of rotate angle

% student ID = 111521035, I = the last number of student ID = 5
beta = mod(5,4)+1;

i = 0:n-1;
alpha = (5*i+beta)*pi/25;

xin = zeros(1,n);
yin = zeros(1,n);
xout = zeros(1,n);
yout = zeros(1,n);
out_angle =zeros(1,n);
phase_error = zeros(n,stage);
float_error = zeros(1,n);

ang = zeros(1,stage);
num = zeros(n,1);

for i = 0:n-1
    xin(1,i+1) = abs(truncation(sin(alpha(1,i+1)), 12));
    yin(1,i+1) = truncation(cos(alpha(1,i+1)), 12);
    out_angle(1,i+1) = atan(yin(1,i+1)/xin(1,i+1));
end

for i =1:n
    num(i,1) = 0;
    cnt = 0;
    mu_i = 0;
    x0 = xin(1,i);
    y0 = yin(1,i);
    ang0 = 0;

%     while(true)
    for k= 1:stage
        mu_i = -sign(y0);

        x1 = x0;
        y1 = y0;
        ang1 = ang0;

        x0 = x1 - truncation(mu_i*(2^(-cnt))*y1, 12);
        y0 = truncation( mu_i*(2^(-cnt))*x1, 12) + y1;
        ang0 = ang1 - mu_i*atan(2^(-cnt));
        
        %
        cnt = cnt+1;
        num(i,1) = num(i,1)+1;

        xout(i,cnt) = x0;
        yout(i,cnt) = y0;
        ang(i,cnt) = ang0;

        % RMSE 中的 Square and Error
        phase_error(i,k) = abs(atan(y0/x0))^2;

        float_error(1,i) = abs(out_angle(1,i)-ang0);

    end
end
% RMSE 中的 Root and Mean
avg_phase_err =sqrt(mean(phase_error));
for avg_cnt = 1:stage
    if((avg_phase_err(1,avg_cnt)) <(0.4*2^(-9)))
        break
    end
end

plot(avg_phase_err)
set(gca, 'YScale', 'log')
title('The average phase error');
xlabel('numbers of micro-rotations 𝑁'), ylabel('the average phase errors')
yline(0.4*2^(-9),'-r','0.4*2 ^-^9')

%% procedure 3_2
% 決定好ang所需的micro-rotation後
% 這裡開始決定 elementary angles的wordlength

clear;clc;

n = 11; % number of input data is 11
stage = 20; % max number of rotate angle
ang_stage = 13;

% student ID = 111521035, I = the last number of student ID = 5
beta = mod(5,4)+1;

i = 0:n-1;
alpha = (5*i+beta)*pi/25;

xin = zeros(1,n);
yin = zeros(1,n);

xout = zeros(1,n);
yout = zeros(1,n);
out_angle =zeros(1,n);
fix_ele_angle = zeros(stage,stage);
float_ele_angle = zeros(stage,stage);

ang = zeros(1,stage);
num = zeros(n,1);

for i = 0:n-1
    % X and Y are quantized into S1.12
    xin(1,i+1) = abs(truncation(sin(alpha(1,i+1)), 12));
    yin(1,i+1) = truncation(cos(alpha(1,i+1)), 12);
    out_angle(1,i+1) = atan(yin(1,i+1)/xin(1,i+1));
end
for k =1:stage
    for i =1:n
        num(i,1) = 0;
        cnt = 0;
        mu_i = 0;
        x0 = xin(1,i);
        y0 = yin(1,i);
        ang0 = 0;
    
          for m = 1:n
            mu_i = -sign(y0);
    
            x1 = x0;
            y1 = y0;
            ang1 = ang0;
    
            x0 = x1 - truncation(mu_i*(2^(-cnt))*y1, 12);
            y0 = truncation( mu_i*(2^(-cnt))*x1, 12) + y1;
            ang0 = ang1 - truncation(mu_i*atan(2^(-cnt)),k);
            %
            cnt = cnt+1;
            num(i,1) = num(i,1)+1;
    
            xout(i,cnt) = x0;
            yout(i,cnt) = y0;
            ang(i,cnt) = ang0;


          end
          % RMSE 中的 SE
          phase_error(i,k) = abs(out_angle(1,i)-ang0)^2;
    end
end
% RMSE 中的 RM
avg_phase_err = sqrt(mean(phase_error));

for avg_cnt = 1:20
    if((avg_phase_err(1,avg_cnt)) <(0.4*2^(-9)))
        break
    end
end

for i = 1:stage         % fraction part
    for j= 1:stage      % N step
        fix_ele_angle(i,j) = truncation(atan(2^(-j+1)),i);
        float_ele_angle(i,j) = atan(2^(-j+1));
    end
end

for i=1:ang_stage+1
    temp = (2^ang_stage)*truncation(fix_ele_angle(ang_stage,i),ang_stage);
    bin_fix = dec2bin(temp);
end

plot(avg_phase_err)
set(gca, 'YScale', 'log')
title('The average phase error versus different elementary angles word-length');
yline(0.4*2^(-9),'-r','0.4*2 ^-^9')

xlabel('bits'), ylabel('the average phase errors')


%% procedure 4
% magnitude error

clear;clc;

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
n = 11; % number of input data is 11
stage = 20; % max number of rotate angle
ang_stage = 13;
% s_n = 0.60725293500924948;

% student ID = 111521035, I = the last number of student ID = 5
beta = mod(5,4)+1;

i = 0:n-1;
alpha = (5*i+beta)*pi/25;

xin = zeros(1,n); 
yin = zeros(1,n);

xout = zeros(n,n); % save x states
yout = zeros(n,n); % save y states
mag_error = zeros(n,stage); % save magnitude errors

ang = zeros(1,stage); % save angle states
num = zeros(n,1); % to count the times that roate the angle

for i = 0:n-1
    xin(1,i+1) = abs(truncation(sin(alpha(1,i+1)), 12));
    yin(1,i+1) = truncation(cos(alpha(1,i+1)), 12);
end

for i =1:n
    
    num(i,1) = 0;
    cnt = 0;
    mu_i = 0;
    x0 = xin(1,i);
    y0 = yin(1,i);
    ang0 = 0;
    origin = sqrt((xin(i)^2)+(yin(i)^2));

    for k= 1:stage
        % to rotate the angle
        mu_i = -sign(y0);

        x1 = x0;
        y1 = y0;
        ang1 = ang0;

        % X and Y are quantized into S1.12
        % x0 = x1 - truncation(mu_i*(2^(-cnt))*y1, 12);
        % y0 = truncation( mu_i*(2^(-cnt))*x1, 12) + y1;
        % ang0 = ang1 - mu_i*atan(2^(-cnt));
        x0 = x1 - (mu_i*(2^(-cnt))*y1);
        y0 = ( mu_i*(2^(-cnt))*x1) + y1;
        ang0 = ang1 - mu_i*atan(2^(-cnt));


        cnt = cnt+1;
        num(i,1) = num(i,1)+1;
        
        % to save the rotation result
        xout(i,cnt) = x0;
        yout(i,cnt) = y0;
        ang(i,cnt) = ang0;
        % RMSE 的 SE
        % mag_error(i,cnt) = ...
        % (origin) - (x0)*s_n(k) / origin;
        mag_error(i,cnt) = ...
        (sqrt( (x0^2)+(y0^2) ) - abs(x0)) / sqrt((x0^2)+(y0^2));
    end
end
% RMSE 的 M
avg_magnitude_err = mean(mag_error);
avg_cnt = 1;
while(true)
    
    if((avg_magnitude_err(1,avg_cnt)) <0.002)
        break
    end
    avg_cnt= avg_cnt+1;
end

plot(avg_magnitude_err)
title('The error of the magnitude versus different number of micro-rotations')
xlabel('the number of the required micro-rotations(N)');
ylabel('error of the magnitude function');
set(gca, 'YScale', 'log')
yline(0.002,'-r','0.2%')

%% procedure 5
% shift-and-add operation

clear;clc;

format long

s_n = ones(1,35);       % s1.12 (14 bits)
temp_2 = zeros(1,35);

for i = 1:35
    temp_2(1,i) = 1/(sqrt(1+2^(-2*(i-1))));
end
for i = 1:35
    for j= 1:i
        s_n(1,i) = s_n(1,i)*(temp_2(1,j));
    end
end

one_of_s_n = 1./s_n;
trun_sn = (2^12) * truncation(s_n, 12);

sn_11 = s_n(1,11); % 使用11個stage的rotation
bin_sn = dec2bin(trun_sn);

the_sn = bin_sn(11,:);
num_sn = zeros(1, length(the_sn)) ;

bin2dec(the_sn)
CSD_in = append('00',the_sn);
CSD_sn = CSD(CSD_in)


N = 35;
temp_2 = zeros(1,N);

n = 11; % number of input data is 11
stage = 20; % max number of rotate angle
ang_stage = 13;

% student ID = 111521035, I = the last number of student ID = 5
beta = mod(5,4)+1;

i = 0:n-1;
alpha = (5*i+beta)*pi/25;

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
    yin(1,i+1) = truncation(cos(alpha(1,i+1)), 12);
    out_angle(1,i+1) = atan(yin(1,i+1)/xin(1,i+1));
end
for k =1:20
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
    
            x0 = x1 - truncation(mu_i*(2^(-cnt))*y1, 12);
            y0 = truncation( mu_i*(2^(-cnt))*x1, 12) + y1;
            ang0 = ang1 - mu_i*atan(2^(-cnt));
            %
            cnt = cnt+1;
            num(i,1) = num(i,1)+1;
    
            xout(i,cnt) = x0;
            yout(i,cnt) = y0;
            ang(i,cnt) = ang0;

          end
          sxout(i,1) = x0 * sn_11;
          sx_compare(i,k) = truncation(x0 * truncation(sn_11,k), k);

          sn_error(i,k) = abs(sxout(i,1)-sx_compare(i,k))^2;
          phase_error(i,k) = abs(out_angle(1,i)-ang0)^2;

    end
end
avg_phase_err = sqrt(mean(sn_error));

for avg_cnt = 1:20
    if((avg_phase_err(1,avg_cnt)) <(0.4*2^(-9)))
        break
    end
end


for i = 1:20 % fraction part
    for j= 1:20 % N step
        fix_ele_angle(i,j) = truncation(atan(2^(-(j-1))),i);
        float_ele_angle(i,j) = atan(2^(-(j-1)));
    end
end

plot(avg_phase_err)
set(gca, 'YScale', 'log')
title('output RMSE versus different S(N) word-length');
yline(0.4*2^(-9),'-r','0.4*2 ^-^9')

xlabel('S(N) word-length(bits)'), ylabel('Output RMSE')
file_ele_angle = fopen('./elementry_ang.txt','w');
file_ele_angle_bin = fopen('./elementry_ang_bin.txt','w');


bin_fix_ele_angle = fix_ele_angle(11,:)*(2^12);

for k =1:11
    nbytes = fprintf(file_ele_angle,'N[%2d] = %d;\n' , k-1, bin_fix_ele_angle(1,k));
    nbytes2 = fprintf(file_ele_angle_bin,'N[%2d] = %s;\n' , k-1, dec2bin(bin_fix_ele_angle(1,k),12));
end

ST = fclose('all');


%% procedure 6

clear
clc


elementary_angles = zeros(1,13);
elementary_angles_out = zeros(1,13);

for i = 1:13
    elementary_angles(1,i) = atan(1/(2^(i-1)));
    elementary_angles_out(1,i) = 2^(13) * truncation( atan(1/(2^(i-1))), 13);
end


%% error between verilog and matlab

error = zeros(1,11);

x1 = linspace(1,11,11);

stem(x1,error)
title('The error between the Verilog output and Matlab output')
xlabel('index of inputs')
ylabel('errors')