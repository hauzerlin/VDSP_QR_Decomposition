function [Q_H, R] = QR_CORDIC_4_4(A,options)
    arguments(Input)
        A 
        options.stage (1,1) {mustBePositive} = 200
        options.input_length (1,1) {mustBeInteger} = 0 %input wordlength
        options.input_angle_length (1,1) {mustBeInteger} = 0
        options.element_angle_length (1,1) {mustBeInteger} = 0
        options.SN_length (1,1) {mustBeInteger} = 0
    end
    SIZE = size(A,1);
    I = eye(SIZE);
    Input_data = cat(2,A,I);
    data = zeros(SIZE,2*SIZE,SIZE);
    data(:,:,1) = Input_data;

    angle_out = zeros(SIZE-1);
    flag_out = zeros(SIZE-1);
% caculate   
    for idk = 1:SIZE-1
        if(idk>1) 
            data(:,:,idk) = data(:,:,idk-1);
        end

        for idx = idk:SIZE-1
            [angle_out(idx-idk+1,idk), data(idk,idk,idk), ...
                data(idx+1,idk,idk), flag_out(idx-idk+1,idk)] = ...
           CORDIC(data(idk,idk,idk), data(idx+1,idk,idk),...
                stage = options.stage, input_length = options.input_length, ...
                input_angle_length = options.input_angle_length , ...
                element_angle_length = options.element_angle_length, ...
                SN_length = options.SN_length);
        end

        for idy = 1+idk:8
            for idx = 1+idk:4
                [~, data(idk,idy,idk), data(idx,idy,idk)] = ...
            CORDIC( data(idk,idy,idk), data(idx,idy,idk),...
                -angle_out(idx-idk,idk), flag_out(idx-idk,idk), 1,...
                stage = options.stage, input_length = options.input_length, ...
                input_angle_length = options.input_angle_length , ...
                element_angle_length = options.element_angle_length, ...
                SN_length = options.SN_length);
            end
        end
        if (idk<SIZE-1)
            for idx = 1:(SIZE-1)-idk
                for idy = 1:(SIZE-idk)-idx
                [~, data(idx+idk,idk,idk), data(idx+idy+idk,idk,idk)] = ...
                    CORDIC(data(idx+idk,idk,idk), data(idx+idy+idk,idk,idk), 0,0,1,...
                stage = options.stage, input_length = options.input_length, ...
                input_angle_length = options.input_angle_length , ...
                element_angle_length = options.element_angle_length, ...
                SN_length = options.SN_length);
                end
            end
        end
    end
    % data(SIZE,:,SIZE) = -data(SIZE,:,SIZE);
    Q_H = data(:,SIZE+1:2*SIZE,SIZE-1);
    R = data(:,1:SIZE,SIZE-1);
end