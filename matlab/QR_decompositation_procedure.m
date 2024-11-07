%% input data
    clc;clear;
    format long
    rng(9487,"twister");

    % set range in (50,100) rand
    % a = 50;
    % b = 100;
    % r = (b-a).*rand(1000,1) + a;
    % A= (1-(-1)) .*rand(4) -1;
    % A = orth(A)';
    % A = [-0.341067858288107 0.128228835541584 -0.338625724761107 -0.241518695807740
    %     -0.780501240379151 0.318329897992542 -0.254694143386287 0.514322409868073
    %     0.0319515900008220 0.00839437248361374 0.836084263956466 0.401688567289137
    %     -0.991013607683554 -0.725621223052855 -0.965140598506631 -0.291424951705322];
    A = [-11926 13869 6827 3702 
         11922 -13810 -7803 13433
         -4665 -6045 -14217 2352
         -14502 11806 -14303 13359];
    A = A/(2^14);
    I = eye(4);
    Input_data = cat(2,A,I);
    Input_data_hard = cat(2,Input_data, zeros(4,3));
    for idx = 2:4
        Input_data_hard(idx,:) = circshift(Input_data_hard(idx,:),idx-1);
    end
    data = zeros(4,8,3);
    data(:,:,1) = Input_data;

% caculate
    % first row
    % caculate first rotate angle for each cloumn.
    
    data(:,:,1) = Input_data;

    for idx = 1:3
        [angle_out(idx,1), data(1,1,1), data(idx+1,1,1), flag_out(idx,1)] = ...
           CORDIC(data(1,1,1), data(idx+1,1,1));
    end

    for idy = 2:8
        for idx = 2:4
            [angle_temp, data(1,idy,1), data(idx,idy,1)] = ...
                CORDIC( data(1,idy,1), data(idx,idy,1), -angle_out(idx-1), flag_out(idx-1), 1);
        end
    end

    for idx = 1:2
        for idy = 1:3-idx
            [~, data(idx+1,1,1), data(idx+idy+1,1,1)] = ...
                CORDIC(data(idx+1,1,1), data(idx+idy+1,1,1), 0,0,1);
        end
    end


    % caculate second ratate angle for each cloumn
    data(:,:,2) = data(:,:,1);

    for idx = 2:3
        [angle_out(idx-1,2), data(2,2,2), data(idx+1,2,2), flag_out(idx-1,2)] = ...
           CORDIC(data(2,2,2), data(idx+1,2,2));
    end

    for idy = 3:8
        for idx = 3:4
            [angle_temp, data(2,idy,2), data(idx,idy,2)] = ...
                CORDIC( data(2,idy,2), data(idx,idy,2), -angle_out(idx-2,2), flag_out(idx-2,2), 1);
        end
    end
    
    for idx = 1:1
        for idy = 1:2-idx
            [~, data(idx+2,2,2), data(idx+idy+2,2,2)] = ...
                CORDIC(data(idx+2,2,2), data(idx+idy+2,2,2), 0,0,1);
        end
    end

    % caculate third ratate angle for each cloumn
    data(:,:,3) = data(:,:,2);

    for idx = 3:3
        [angle_out(idx-2,3), data(3,3,3), data(idx+1,3,3), flag_out(idx-2,3)] = ...
           CORDIC(data(3,3,3), data(idx+1,3,3));
    end

    for idy = 4:8
        for idx = 4:4
            [angle_temp, data(3,idy,3), data(idx,idy,3)] = ...
                CORDIC( data(3,idy,3), data(idx,idy,3), -angle_out(idx-3,3), flag_out(idx-3,3), 1);
        end
    end    
    % caculate layer 4
    data(4,:,3) = -data(4,:,3);

% error caculate

    Q_H = data(:,5:8,3);
    Q = Q_H';
    R = data(:,1:4,3);
    
    A_QR = Q*R;

    QR_error = (A_QR-A);

    % A_binary = truncation(A,input_length) .* (2^14)
    % I_binary = eye(4) .* (2^14)
    % 
    % Q_binary = Q_H * (2^14)
    % R_binary = R * (2^14)

%% Test function QR_CORDIC_4_4
    clc;clear;
    format long
    rng(9487,"twister");
    stage = 13;
    input_length = 14; input_angle_length =14;
    element_angle_length = 14; SN_length = 14;
    % A= (1-(-1)) .*rand(4) -1;
    % A = orth(A)';
    A = [-11926 13869 6827 3702 
     11922 -13810 -7803 13433
     -4665 -6045 -14217 2352
     -14502 11806 -14303 13359];
    A = A/(2^14);
    [Q_test, R_test] = QR_CORDIC_4_4(A,...
                stage = stage, input_length = input_length, ...
                input_angle_length = input_angle_length , ...
                element_angle_length = element_angle_length, ...
                SN_length = SN_length);
    A_QR_test = Q_test'*R_test;

    QR_error_test = (A_QR_test-A);
    QR_RMSE_test = norm(A_QR_test-A)/norm(A);
    
    A_binary = truncation(A,input_length) .* (2^14)
    I_binary = eye(4) .* (2^14)

    Q_binary = Q_test * (2^14)
    R_binary = R_test * (2^14)

    % for verilog pattern
    pattern1 = cat(1,A_binary(:,1),I_binary(:,1));
    pattern2 = cat(1,A_binary(:,2),I_binary(:,2));
    pattern3 = cat(1,A_binary(:,3),I_binary(:,3));
    pattern4 = cat(1,A_binary(:,4),I_binary(:,4));

    in1 = fopen("hower\in1.txt","w");
    in2 = fopen("hower\in2.txt","w");
    in3 = fopen("hower\in3.txt","w");
    in4 = fopen("hower\in4.txt","w");

    for idx = 1:length(pattern1)
        fprintf(in1,"%s\n", dec2bin(pattern1(idx),input_length+2));
        fprintf(in2,"%s\n", dec2bin(pattern2(idx),input_length+2));
        fprintf(in3,"%s\n", dec2bin(pattern3(idx),input_length+2));
        fprintf(in4,"%s\n", dec2bin(pattern4(idx),input_length+2));
    end
    rng(9453,"twister");
    stage = 13;
    input_length = 14; input_angle_length =14;
    element_angle_length = 14; SN_length = 14;
    A= (1-(-1)) .*rand(4) -1;
    A = orth(A)';
    [Q_test, R_test] = QR_CORDIC_4_4(A,...
                stage = stage, input_length = input_length, ...
                input_angle_length = input_angle_length , ...
                element_angle_length = element_angle_length, ...
                SN_length = SN_length);
    A_QR_test = Q_test'*R_test;

    QR_error_test = (A_QR_test-A);
    QR_RMSE_test = norm(A_QR_test-A)/norm(A);

    A_binary = truncation(A,input_length) .* (2^14)
    I_binary = eye(4) .* (2^14)

    Q_binary = Q_test * (2^14)
    R_binary = R_test * (2^14)

    % for verilog pattern
    pattern1 = cat(1,A_binary(:,1),I_binary(:,1));
    pattern2 = cat(1,A_binary(:,2),I_binary(:,2));
    pattern3 = cat(1,A_binary(:,3),I_binary(:,3));
    pattern4 = cat(1,A_binary(:,4),I_binary(:,4));

    for idx = 1:length(pattern1)
        fprintf(in1,"%s\n", dec2bin(pattern1(idx),input_length+2));
        fprintf(in2,"%s\n", dec2bin(pattern2(idx),input_length+2));
        fprintf(in3,"%s\n", dec2bin(pattern3(idx),input_length+2));
        fprintf(in4,"%s\n", dec2bin(pattern4(idx),input_length+2));
    end
        rng(9981,"twister");
    stage = 13;
    input_length = 14; input_angle_length =14;
    element_angle_length = 14; SN_length = 14;
    A= (1-(-1)) .*rand(4) -1;
    A = orth(A)';
    [Q_test, R_test] = QR_CORDIC_4_4(A,...
                stage = stage, input_length = input_length, ...
                input_angle_length = input_angle_length , ...
                element_angle_length = element_angle_length, ...
                SN_length = SN_length);
    A_QR_test = Q_test'*R_test;

    QR_error_test = (A_QR_test-A);
    QR_RMSE_test = norm(A_QR_test-A)/norm(A);

    A_binary = truncation(A,input_length) .* (2^14)
    I_binary = eye(4) .* (2^14)

    Q_binary = Q_test * (2^14)
    R_binary = R_test * (2^14)

    % for verilog pattern
    pattern1 = cat(1,A_binary(:,1),I_binary(:,1));
    pattern2 = cat(1,A_binary(:,2),I_binary(:,2));
    pattern3 = cat(1,A_binary(:,3),I_binary(:,3));
    pattern4 = cat(1,A_binary(:,4),I_binary(:,4));

    for idx = 1:length(pattern1)
        fprintf(in1,"%s\n", dec2bin(pattern1(idx),input_length+2));
        fprintf(in2,"%s\n", dec2bin(pattern2(idx),input_length+2));
        fprintf(in3,"%s\n", dec2bin(pattern3(idx),input_length+2));
        fprintf(in4,"%s\n", dec2bin(pattern4(idx),input_length+2));
    end

    fclose('all');
%% gerenate pattern for quntaziation

    clc;clear;
    format long
    rng(9487,"twister");
    A= (1-(-1)) .*rand(4) -1;
    A = orth(A);
    rotation_number = 50;
tic
for idx = 1:rotation_number
    [Q_test(:,:,idx), R_test(:,:,idx)] = QR_CORDIC_4_4(A, stage=idx);
    A_QR_test = Q_test(:,:,idx)'*R_test(:,:,idx);

    QR_error_test(:,:,idx) = (A_QR_test-A);
    QR_RMSE_test(idx) = sqrt(mean(mean((A_QR_test-A).^2)'));
end
toc

plot(1:rotation_number,20*log10(QR_RMSE_test));
% function [Q_H, R] = QR_CORDIC_4_4(A,options)
%     arguments(Input)
%         A (4,4)
%         options.stage (1,1) {mustBePositive} = 200
%         options.input_length (1,1) {mustBeInteger} = 0 %input wordlength
%         options.input_angle_length (1,1) {mustBeInteger} = 0
%         options.element_angle_length (1,1) {mustBeInteger} = 0
%         options.SN_length (1,1) {mustBeInteger} = 0
%     end
%     SIZE = size(A,1);
%     I = eye(4);
%     Input_data = cat(2,A,I);
%     data = zeros(SIZE,2*SIZE,SIZE);
%     data(:,:,1) = Input_data;
% 
%     angle_out = zeros(SIZE-1);
%     flag_out = zeros(SIZE-1);
% % caculate   
%     for idk = 1:SIZE-1
%         if(idk>1) 
%             data(:,:,idk) = data(:,:,idk-1);
%         end
% 
%         for idx = idk:SIZE-1
%             [angle_out(idx-idk+1,idk), data(idk,idk,idk), ...
%                 data(idx+1,idk,idk), flag_out(idx-idk+1,idk)] = ...
%            CORDIC(data(idk,idk,idk), data(idx+1,idk,idk),...
%                 stage = options.stage, input_length = options.input_length, ...
%                 input_angle_length = options.input_angle_length , ...
%                 element_angle_length = options.element_angle_length, ...
%                 SN_length = options.SN_length);
%         end
% 
%         for idy = 1+idk:8
%             for idx = 1+idk:4
%                 [~, data(idk,idy,idk), data(idx,idy,idk)] = ...
%             CORDIC( data(idk,idy,idk), data(idx,idy,idk),...
%                 -angle_out(idx-idk,idk), flag_out(idx-idk,idk), 1,...
%                 stage = options.stage, input_length = options.input_length, ...
%                 input_angle_length = options.input_angle_length , ...
%                 element_angle_length = options.element_angle_length, ...
%                 SN_length = options.SN_length);
%             end
%         end
%     end
%     Q_H = data(:,SIZE+1:2*SIZE,SIZE-1);
%     R = data(:,1:SIZE,SIZE-1);
% end