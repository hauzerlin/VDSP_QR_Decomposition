%% Canonic Signed Digit (CSD) function
function Ci = CSD_wayne(B)

% B = [0 0 0 1 0 1 1 0];
% B = [0 0 1 1 1 0 1 0 1 1 1];
% B = [1 1 0 0 1 0 1 1 1 0 1];
% 以上是測試用的值

Bi = [1 B 0];       % 初始化 : 設定 Bi 頭尾的值

Si = [];
Gi = [];

L = length(Bi);

Gi(L) = 0;          % 初始化 : 設定 Gi 最後面的值

for i = 2 : L-1
    Si(i) = xor(Bi(i), Bi(i+1));
end

for i = 2 : L-1
    if(Gi(L-i+2)>0)
        Gi(L-i+1) = 0;
    else
        Gi(L-i+1) = Si(L-i+1);
    end
end

Ci = [];

for i = 1 : L-2
    Ci(i) = (1-2*Bi(i))*Gi(i+1);
end

end
