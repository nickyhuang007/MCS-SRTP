function[data,arch,cal_mean,cal_diff,cal_SSE] = cal_garch(Idx,constant,Garch,Arch,x0)


% Idx = xlsread('C:\Users\Lenovo\Desktop\国创\new_000001.SS.xlsx','E7001 : E7238');
Idx = Idx';

new_Idx = zeros(length(Idx),1);
for i = 2 : length(Idx)
    new_Idx(i) = (Idx(i) - Idx(i-1))/Idx(i-1);
end

h0 = x0 / normrnd(0,1);

data = zeros(1,length(Idx));
arch = zeros(1,length(Idx));

for i = 1 : length(Idx)
    h1 = sqrt(constant + Garch * x0^2 + Arch * h0^2);
    arch(1,i) = h1^2;
    x1 = h1 * normrnd(0,1);
    data(1,i) = 0 + x1;
    x0 = x1;
    h0 = h1;
end

% [lbq_data_result,lbq_data_p,lbq_data,lbq_data_accept] = lbqtest(data,'Lags',20); %lbq测试,原假设为与后面20位无关,彼此独立 

compare = zeros(1,length(Idx));
new_Idx = new_Idx';
for i = 1 : length(Idx)
    compare(1,i) = new_Idx(1,i) - data(1,i);
end

count = 0;
for i = 1 : length(Idx)
    if(abs(compare(1,i)) < 1)
        count = count + 1;
    end
end

x = 1 : length(Idx);
cal_mean = mean(new_Idx);
tem_cal_diff = 0;
for i = 1 : length(new_Idx)
   tem_cal_diff = tem_cal_diff + (new_Idx(i) - cal_mean)^2;
end
tem_cal_diff = tem_cal_diff/length(new_Idx);
tem_cal_SSE = 0;
for i = 1 : length(new_Idx)
    tem_cal_SSE = tem_cal_SSE + compare(i)^2;
end

cal_diff = tem_cal_diff;
cal_SSE = tem_cal_SSE;