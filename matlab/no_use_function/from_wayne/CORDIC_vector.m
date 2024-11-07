%% Generate the test pattern
clc
clear

% Parameter
ROTATION_NUM_FIX = 10;
WORDLENGTH_PHASE = 12;
WORDLENGTH_MAG = 12;
WORDLENGTH_SF = 12;
PATTERN_NUM = 10000;
ID_LAST_DIGIT = 0;

% alpha
beta = mod(ID_LAST_DIGIT, 4) + 1;
n = linspace(0, PATTERN_NUM-1 ,PATTERN_NUM);
alpha = 2*pi*rand(1, PATTERN_NUM);
%%
% X and Y (Cartesian coordinate)
X = sin(alpha);
Y = cos(alpha);

X_fixed12 = truncate(X, WORDLENGTH_MAG);
Y_fixed12 = truncate(Y, WORDLENGTH_MAG);

% Initial stage (Map to 1/4 quadrants)
X_map = zeros(1, PATTERN_NUM);
Y_map = zeros(1, PATTERN_NUM);
quadrant = zeros(1, PATTERN_NUM);

for idx = 1:PATTERN_NUM
    if(X_fixed12(idx)>=0)
        if(Y_fixed12(idx)>=0)
            X_map(idx) = X_fixed12(idx);
            Y_map(idx) = Y_fixed12(idx);
            quadrant(idx) = 1;
        else
            X_map(idx) = X_fixed12(idx);
            Y_map(idx) = Y_fixed12(idx);
            quadrant(idx) = 4;
        end
    elseif (X_fixed12(idx)<0)
        if(Y_fixed12(idx)>=0)
            X_map(idx) = Y_fixed12(idx);
            Y_map(idx) = -X_fixed12(idx);
            quadrant(idx) = 2;
        else
            X_map(idx) = -Y_fixed12(idx);
            Y_map(idx) = X_fixed12(idx);
            quadrant(idx) = 3;
        end
    end
end
%%
% Calculate X(N) and Y(N), 0<N<11 by CORDIC operation (fixed point)
X_rot_fixed_12 = zeros(ROTATION_NUM_FIX + 1, PATTERN_NUM);
Y_rot_fixed_12 = zeros(ROTATION_NUM_FIX + 1, PATTERN_NUM);
X_rot10_fixed_12 = zeros(1, PATTERN_NUM);
Y_rot10_fixed_12 = zeros(1, PATTERN_NUM);

theta_rot_fixed_12 = zeros(ROTATION_NUM_FIX + 1, PATTERN_NUM);
theta_rot10_fixed_12 = zeros(1, PATTERN_NUM);

% Elemantary angle(S1.12)
ele_angle_fixed_12 = zeros(1, ROTATION_NUM_FIX);

for idx = 1:ROTATION_NUM_FIX
    e1 = atan(2^(-(idx-1)));
    ele_angle_fixed_12(idx) = truncate(e1, WORDLENGTH_PHASE);
end

% scaling factor 𝑆(𝑁)
SF = zeros(1, ROTATION_NUM_FIX);
% Calculate scaling factor
for idx = 1:ROTATION_NUM_FIX
    if (idx == 1)
        SF(idx) = 1/sqrt(1+2^(-2*(idx-1)));
    else
        SF(idx) = SF(idx-1) * 1/sqrt(1+2^(-2*(idx-1)));
    end
end

% CORDIC operation (fixed point)
for idx = 1:ROTATION_NUM_FIX + 1
    if (idx == 1)                   % initialization
        X_rot_fixed_12(idx, :) = truncate(X_map, WORDLENGTH_MAG);
        Y_rot_fixed_12(idx, :) = truncate(Y_map, WORDLENGTH_MAG);
        theta_rot_fixed_12(idx, :) = 0;
    else
        for idy = 1:PATTERN_NUM
            if(Y_rot_fixed_12(idx-1, idy) >= 0)
                mu(idy) = -1;
            else
                mu(idy) = 1;
            end
        end
%         b = X_rot_fixed_12(idx-1, :) - mu.*2^(-(idx-2)).*Y_rot_fixed_12(idx-1, :);
        b1 = 2^(-(idx-2));
        b1_t = truncate(b1, WORDLENGTH_MAG);
        b2 = b1_t.*Y_rot_fixed_12(idx-1, :);
        b2_t = truncate(b2, WORDLENGTH_MAG);
        b3 = X_rot_fixed_12(idx-1, :) - (mu.*b2_t);
        X_rot_fixed_12(idx, :) = truncate(b3, WORDLENGTH_MAG);
        
%         c = mu.*2^(-(idx-2)).*X_rot_fixed_12(idx-1, :) + Y_rot_fixed_12(idx-1, :);
        c1 = 2^(-(idx-2));
        c1_t = truncate(c1, WORDLENGTH_MAG);
        c2 = c1_t.*X_rot_fixed_12(idx-1, :);
        c2_t = truncate(c2, WORDLENGTH_MAG);
        c3 = (mu.*c2_t) + Y_rot_fixed_12(idx-1, :);
        Y_rot_fixed_12(idx, :) = truncate(c3, WORDLENGTH_MAG);
        
%         d = theta_rot_fixed_12 - mu.*ele_angle_fixed_12(idx-1);
        d1 = mu.*ele_angle_fixed_12(idx-1);
        d1_t = truncate(d1, WORDLENGTH_PHASE);
        d2 = theta_rot_fixed_12(idx-1, :) - d1_t;
        theta_rot_fixed_12(idx, :) = truncate(d2, WORDLENGTH_PHASE);
    end
end

% X(10), Y(10), 𝜃(10) by fixed point
X_rot10_fixed_12 = X_rot_fixed_12(ROTATION_NUM_FIX + 1, :);
Y_rot10_fixed_12 = Y_rot_fixed_12(ROTATION_NUM_FIX + 1, :);
theta_rot10_fixed_12 = theta_rot_fixed_12(ROTATION_NUM_FIX + 1, :);

% Scaling
% g = SF_10_fixed_12 * X_rot10_fixed_12;
g1 = truncate(2^(-1), WORDLENGTH_SF);
g2 = truncate(2^(-3), WORDLENGTH_SF);
g3 = truncate(2^(-6), WORDLENGTH_SF);
g4 = truncate(2^(-9), WORDLENGTH_SF);
g5 = truncate(2^(-12), WORDLENGTH_SF);

% h = X_rot10_fixed_12 * (g1 + g2 - g3 - g4 - g5);
h1 = truncate(X_rot10_fixed_12 * g1, WORDLENGTH_SF);
h2 = truncate(X_rot10_fixed_12 * g2, WORDLENGTH_SF);
h3 = truncate(X_rot10_fixed_12 * g3, WORDLENGTH_SF);
h4 = truncate(X_rot10_fixed_12 * g4, WORDLENGTH_SF);
h5 = truncate(X_rot10_fixed_12 * g5, WORDLENGTH_SF);

m1 = truncate(h1 + h2, WORDLENGTH_SF);
m2 = truncate(h3 + h4, WORDLENGTH_SF);
m3 = truncate(m1 - m2, WORDLENGTH_SF);
m4 = truncate(m3 - h5, WORDLENGTH_SF);
X_rot10_fixed12_scaling = truncate(m4, WORDLENGTH_SF);