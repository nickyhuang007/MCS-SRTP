clc
clear

% Idx = xlsread('C:\Users\Lenovo\Desktop\国创\道琼斯工业平均指数.xls','E2 : E22754');%第一位为0，因为为第一天手动置为0
% [~,Trd] = xlsread('C:\Users\Lenovo\Desktop\国创\道琼斯工业平均指数.xls','C2 : C22754');
Idx = xlsread('C:\Users\Lenovo\Desktop\国创\new_000001.SS.xlsx','E2 : E7000');
[~,Trd] = xlsread('C:\Users\Lenovo\Desktop\国创\new_000001.SS.xlsx','A2 : A7000');

Idx = Idx';

[year,month,day] = get_date(Trd);

data = zeros(length(Idx),4);
data(:,1) = year(:,1);
data(:,2) = month(:,1);
data(:,3) = day(:,1);
data(:,4) = Idx(1,:);

Idx_year = zeros(1,100);
count_year = zeros(1,100);

for i = 1 : length(Idx)
    Idx_year(1,data(i,1) - 1928 + 1) = Idx_year(1,data(i,1) - 1928 + 1) + data(i,4);
    count_year(1,data(i,1) - 1928 + 1) = count_year(1,data(i,1) - 1928 + 1) + 1;
end

i = 1;
while(count_year(1,i) ~= 0)    
    Idx_year(1,i) = Idx_year(1,i) / count_year(1,i);
    i = i + 1;
end

new_Idx = zeros(length(Idx),1);
for i = 2 : length(Idx)
    new_Idx(i) = (Idx(i) - Idx(i-1))/Idx(i-1);
end

average_Idx = 0;
for i = 2 : length(Idx)
    average_Idx = average_Idx + Idx(1,i);
end
average_Idx = average_Idx/length(Idx); %均值

standard_Idx = std(Idx,0,2); %标准差

skew_Idx = skewness(Idx(1,:)); %偏度

kurt_Idx = kurtosis(Idx(1,:)); %峰度

[lbq_Idx_result,lbq_Idx_p,lbq_Idx,lbq_Idx_accept] = lbqtest(Idx,'Lags',20); %lbq测试,原假设为与后面20位无关,彼此独立

[arch_Idx_result,arch_Idx_p,arch_Idx,arch_Idx_accept] = archtest(Idx);

write_data = [average_Idx standard_Idx skew_Idx kurt_Idx lbq_Idx arch_Idx];

xlswrite('C:\Users\Lenovo\Desktop\国创\收益率基本统计特性.xlsx',write_data,'B2 : G2');

fun = garch(1,1);

funfit = estimate(fun,new_Idx);

fun2 = gjr(1,1);

funfit2 = estimate(fun2,new_Idx);

fun3 = egarch(1,1);

funfit3 = estimate(fun3,new_Idx);
% write_data = [funfit.Offset funfit.Constant funfit.ARCH{1,1}(1,1) funfit.GARCH{1,1}(1,1)];
% xlswrite('C:\Users\Lenovo\Desktop\国创\Garch(1,1)参数.xlsx',write_data,'E11 : H11');
