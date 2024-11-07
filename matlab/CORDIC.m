function [angle_out, x_out, y_out, flag_out] = ...
         CORDIC(x, y, angle_in, flag_in, mode, options)
% mode 0:vectoring mode(clockwise), mode 1:rotation mode(anticlockwise)
% offical use for final project

    arguments(Input)
        x (1,:) double
        y (1,:) double
        angle_in double = 0
        flag_in {mustBeInteger} = 0
        mode (1,1) = 0 
        options.stage (1,1) {mustBePositive} = 200
        options.input_length (1,1) {mustBeInteger} = 0 %input wordlength
        options.input_angle_length (1,1) {mustBeInteger} = 0
        options.element_angle_length (1,1) {mustBeInteger} = 0
        options.SN_length (1,1) {mustBeInteger} = 0
    end
    arguments(Output)
        angle_out
        x_out
        y_out
        flag_out
    end
    format long

    % S(N) = 1/Π(1+2^(-2i))
    s_n = ones(1,options.stage+1);
    temp = zeros(1,options.stage+1);
    for i = 1:options.stage+1
        temp(1,i) = 1/(sqrt(1+2^(-2*(i-1)))); 
    end
    for i = 1:options.stage+1
        for j= 1:i
            s_n(1,i) = s_n(1,i)*(temp(1,j));
        end
    end
    s_n = truncation(s_n,options.SN_length);

    % elementary angle 
    elementary_angles = zeros(1,options.stage+1);
    for i = 1:options.stage+1
        elementary_angles(1,i) =  ...
            truncation(atan(1/(2^(i-1))), options.element_angle_length);
    end

    % initial options.stage
    x_data = zeros(1,options.stage+1);
    y_data = zeros(1,options.stage+1);
    angle_data = zeros(1,options.stage+1);
    init_angle = truncation(angle_in,options.input_angle_length);
    
        if(mode==0)
            if(x<0)
                mapping_flag = 1;
            else
                mapping_flag = 0;
            end
            x_data(1) = abs(truncation(x,options.input_length));
            y_data(1) = truncation(y,options.input_length);
            init_angle = 0;
            angle_data(1) = init_angle;
    
        elseif (mode==1)
            if(flag_in==1)
                x_data(1) = -truncation(x,options.input_length);
                y_data(1) = -truncation(y,options.input_length);
            else
                x_data(1) = truncation(x,options.input_length);
                y_data(1) = truncation(y,options.input_length);
            end
            angle_data(1) = init_angle;
            mapping_flag = flag_in;
        else 
            error("mode is wrong");
        end
    % x_data(1) * (2^options.input_length)
    % y_data(1) * (2^options.input_length)
    % caculate stage
    for idx = 1:options.stage
        if (mode==0) 
            if(y_data(idx)<0) mu = 1;
            else mu = -1;
            end
            % mu = -sign(y_data(idx));
        else 
            if(angle_data(idx)<0) mu=-1;
            else mu = 1;
            end
            % mu = sign(angle_data(idx));
        end
        % mu
        x_data(idx+1) = truncation (x_data(idx) - ...
                        mu*...
                        truncation( ...
                        ( truncation( (2^(-(idx-1)) ),options.input_length)* ...
                        y_data(idx) ),options.input_length) ,options.input_length);
        y_data(idx+1) = truncation(mu* ...
                        truncation( (x_data(idx)* ...
                        truncation( (2^(-(idx-1)) ),options.input_length) ), ...
                        options.input_length ) ...
                        + y_data(idx) ,options.input_length);
        if(x_data(idx)~=0 || y_data(idx)~=0)
            angle_data(idx+1) = angle_data(idx) - mu*elementary_angles(idx);
        else
            angle_data(idx+1) = 0;
        end
    end
    % 
    % x_data(2) * (2^options.input_length)
    % y_data(2) * (2^options.input_length)
    x_data(options.stage+1) * (2^options.input_length)
    y_data(options.stage+1) * (2^options.input_length)

    % final stage(output): CSD and mapping
    if(mode==0) % vectroing mode output
        if(mapping_flag==1)
            angle_out = -angle_data(options.stage+1);
        else
            angle_out = angle_data(options.stage+1);
        end
        if(options.SN_length==0)
            x_out = x_data(options.stage+1) * s_n(options.stage);
            y_out = y_data(options.stage+1) * s_n(options.stage);
        else
            CSD_data = CSD( ['0', dec2bin((2^options.SN_length)* ...
                truncation(s_n(options.stage),options.SN_length))] );
            x_CSD = 0; y_CSD = 0;
            for idx =2:length(CSD_data)
                if(CSD_data(idx)==1)
                   x_CSD = x_CSD + ...
                       truncation(x_data(options.stage+1)* ...
                       truncation(2^-(idx-1),options.input_length),options.input_length); 
                   y_CSD = y_CSD + ...
                       truncation(y_data(options.stage+1)* ...
                       truncation(2^-(idx-1),options.input_length),options.input_length); 
                elseif(CSD_data(idx)==-1)
                   x_CSD = x_CSD - ...
                       truncation(x_data(options.stage+1)* ...
                       truncation(2^-(idx-1),options.input_length),options.input_length); 
                   y_CSD = y_CSD - ...
                       truncation(y_data(options.stage+1)* ...
                       truncation(2^-(idx-1),options.input_length),options.input_length); 
                end
            end
            x_out = truncation(x_CSD,options.input_length);
            y_out = truncation(y_CSD,options.input_length);
        end
        flag_out  = mapping_flag;
    
    elseif (mode==1) % rotation mode output
        angle_out = angle_data(options.stage+1);

        % x_out = x_data(options.stage+1) * s_n(options.stage);
        % y_out = y_data(options.stage+1) * s_n(options.stage);
        if(options.SN_length==0)
            x_out = x_data(options.stage+1) * s_n(options.stage);
            y_out = y_data(options.stage+1) * s_n(options.stage);
        else
            CSD_data = CSD( ['0', dec2bin((2^options.SN_length)* ...
                truncation(s_n(options.stage),options.SN_length))] );
            x_CSD = 0; y_CSD = 0;
            for idx =2:length(CSD_data)
                if(CSD_data(idx)==1)
                   x_CSD = x_CSD + ...
                       truncation(x_data(options.stage+1)* ...
                       truncation(2^-(idx-1),options.input_length),options.input_length); 
                   y_CSD = y_CSD + ...
                       truncation(y_data(options.stage+1)* ...
                       truncation(2^-(idx-1),options.input_length),options.input_length); 
                elseif(CSD_data(idx)==-1)
                   x_CSD = x_CSD - ...
                       truncation(x_data(options.stage+1)* ...
                       truncation(2^-(idx-1),options.input_length),options.input_length); 
                   y_CSD = y_CSD - ...
                       truncation(y_data(options.stage+1)* ...
                       truncation(2^-(idx-1),options.input_length),options.input_length); 
                end
            end
            x_out = truncation(x_CSD,options.input_length);
            y_out = truncation(y_CSD,options.input_length);
        end

        flag_out  = flag_in;
    end
% s_n(options.stage)
end