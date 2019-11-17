function[Tr,Tsq] = get_statistics(d_MSE,d_MAE,d_HMSE,d_HMAE,d_QLIKE,d_R2LOG,M0,M)

%输入损失函数的差size = (M0,M0,M)，并输出统计量Tr、Tsq

%d_MSE,d_MAE,d_HMSE,d_HMAE,d_QLIKE,d_R2LOG记录3个模型分别于另外两个模型之间的差,分6种函数

d_mean_MSE = zeros(M0,M0);
d_mean_MAE = zeros(M0,M0);
d_mean_HMSE = zeros(M0,M0);
d_mean_HMAE = zeros(M0,M0);
d_mean_QLIKE = zeros(M0,M0);
d_mean_R2LOG = zeros(M0,M0);

for i = 1 : M0
    for j = 1 : M0
        for k = 1 : M
            d_mean_MSE(i,j) = d_mean_MSE(i,j) + d_MSE(i,j,k);
            d_mean_MAE(i,j) = d_mean_MAE(i,j) + d_MAE(i,j,k);
            d_mean_HMSE(i,j) = d_mean_HMSE(i,j) + d_HMSE(i,j,k);
            d_mean_HMAE(i,j) = d_mean_HMAE(i,j) + d_HMAE(i,j,k);
            d_mean_QLIKE(i,j) = d_mean_QLIKE(i,j) + d_QLIKE(i,j,k);
            d_mean_R2LOG(i,j) = d_mean_R2LOG(i,j) + d_R2LOG(i,j,k);
        end
    end
end

d_mean_MSE = d_mean_MSE/M;
d_mean_MAE = d_mean_MAE/M;
d_mean_HMSE = d_mean_HMSE/M;
d_mean_HMAE = d_mean_HMAE/M;
d_mean_QLIKE = d_mean_QLIKE/M;
d_mean_R2LOG = d_mean_R2LOG/M;

d_var_MSE = zeros(M0,M0); 
d_var_MAE = zeros(M0,M0);
d_var_HMSE = zeros(M0,M0);
d_var_HMAE = zeros(M0,M0);
d_var_QLIKE = zeros(M0,M0);
d_var_R2LOG = zeros(M0,M0);

for i = 1 : M0
    for j = 1 : M0
        for k = 1 : M
            d_var_MSE(i,j) = d_var_MSE(i,j) + (d_MSE(i,j,k) - d_mean_MSE(i,j))^2;
            d_var_MAE(i,j) = d_var_MAE(i,j) + (d_MAE(i,j,k) - d_mean_MAE(i,j))^2;
            d_var_HMSE(i,j) = d_var_HMSE(i,j) + (d_HMSE(i,j,k) - d_mean_HMSE(i,j))^2;
            d_var_HMAE(i,j) = d_var_HMAE(i,j) + (d_HMAE(i,j,k) - d_mean_HMAE(i,j))^2;
            d_var_QLIKE(i,j) = d_var_QLIKE(i,j) + (d_QLIKE(i,j,k) - d_mean_QLIKE(i,j))^2;
            d_var_R2LOG(i,j) = d_var_R2LOG(i,j) + (d_R2LOG(i,j,k) - d_mean_R2LOG(i,j))^2;
        end
    end
end

d_var_MSE = d_var_MSE/M;
d_var_MAE = d_var_MAE/M;
d_var_HMSE = d_var_HMSE/M;
d_var_HMAE = d_var_HMAE/M;
d_var_QLIKE = d_var_QLIKE/M;
d_var_R2LOG = d_var_R2LOG/M;







Tr = zeros(M0,6);
Tsq = zeros(M0,6);

temp_Tr = zeros(M0,M0,6);
temp_Tsq = zeros(M0,M0,6);

for i = 1 : M0
    for j = 1 : M0
        temp_Tr(i,j,1) = (abs(d_mean_MSE(i,j))/sqrt(d_var_MSE(i,j)));
        temp_Tr(i,j,2) = (abs(d_mean_MAE(i,j))/sqrt(d_var_MAE(i,j)));
        temp_Tr(i,j,3) = (abs(d_mean_HMSE(i,j))/sqrt(d_var_HMSE(i,j)));
        temp_Tr(i,j,4) = (abs(d_mean_HMAE(i,j))/sqrt(d_var_HMAE(i,j)));
        temp_Tr(i,j,5) = (abs(d_mean_QLIKE(i,j))/sqrt(d_var_QLIKE(i,j)));
        temp_Tr(i,j,6) = (abs(d_mean_R2LOG(i,j))/sqrt(d_var_R2LOG(i,j)));
    end
end

for i = 1 : M0
    Tr(i,1) = max(temp_Tr(i,:,1));
    Tr(i,2) = max(temp_Tr(i,:,2));
    Tr(i,3) = max(temp_Tr(i,:,3));
    Tr(i,4) = max(temp_Tr(i,:,4));
    Tr(i,5) = max(temp_Tr(i,:,5));
    Tr(i,6) = max(temp_Tr(i,:,6));
end

d_mean_mean_MSE = zeros(1,M0);
d_mean_mean_MAE = zeros(1,M0);
d_mean_mean_HMSE = zeros(1,M0);
d_mean_mean_HMAE = zeros(1,M0);
d_mean_mean_QLIKE = zeros(1,M0);
d_mean_mean_R2LOG = zeros(1,M0);

for i = 1 : M0
    for j = 1 : M0
        if(d_mean_MSE(i,j) ~= 0)
           d_mean_mean_MSE(1,i) = d_mean_mean_MSE(1,i) + d_mean_MSE(i,j);
           d_mean_mean_MAE(1,i) = d_mean_mean_MAE(1,i) + d_mean_MAE(i,j);
           d_mean_mean_HMSE(1,i) = d_mean_mean_HMSE(1,i) + d_mean_HMSE(i,j);
           d_mean_mean_HMAE(1,i) = d_mean_mean_HMAE(1,i) + d_mean_HMAE(i,j);
           d_mean_mean_QLIKE(1,i) = d_mean_mean_QLIKE(1,i) + d_mean_QLIKE(i,j);
           d_mean_mean_R2LOG(1,i) = d_mean_mean_R2LOG(1,i) + d_mean_R2LOG(i,j);
        end
    end
end

d_mean_mean_MSE = d_mean_mean_MSE/(M0-1);
d_mean_mean_MAE = d_mean_mean_MAE/(M0-1);
d_mean_mean_HMSE = d_mean_mean_HMSE/(M0-1);
d_mean_mean_HMAE = d_mean_mean_HMAE/(M0-1);
d_mean_mean_QLIKE = d_mean_mean_QLIKE/(M0-1);
d_mean_mean_R2LOG = d_mean_mean_R2LOG/(M0-1);

d_mean_var_MSE = zeros(1,M0);
d_mean_var_MAE = zeros(1,M0);
d_mean_var_HMSE = zeros(1,M0);
d_mean_var_HMAE = zeros(1,M0);
d_mean_var_QLIKE = zeros(1,M0);
d_mean_var_R2LOG = zeros(1,M0);

for i = 1 : M0
    for j = 1 : M0
        if(d_mean_MSE(i,j) ~= 0)
           d_mean_var_MSE(1,i) = d_mean_var_MSE(1,i) + (d_mean_MSE(i,j) - d_mean_mean_MSE(1,i))^2;
           d_mean_var_MAE(1,i) = d_mean_var_MAE(1,i) + (d_mean_MAE(i,j) - d_mean_mean_MAE(1,i))^2;
           d_mean_var_HMSE(1,i) = d_mean_var_HMSE(1,i) + (d_mean_HMSE(i,j) - d_mean_mean_HMSE(1,i))^2;
           d_mean_var_HMAE(1,i) = d_mean_var_HMAE(1,i) + (d_mean_HMAE(i,j) - d_mean_mean_HMAE(1,i))^2;
           d_mean_var_QLIKE(1,i) = d_mean_var_QLIKE(1,i) + (d_mean_QLIKE(i,j) - d_mean_mean_QLIKE(1,i))^2;
           d_mean_var_R2LOG(1,i) = d_mean_var_R2LOG(1,i) + (d_mean_R2LOG(i,j) - d_mean_mean_R2LOG(1,i))^2;
        end
    end
end

d_mean_var_MSE = d_mean_var_MSE/(M0-1);
d_mean_var_MAE = d_mean_var_MAE/(M0-1);
d_mean_var_HMSE = d_mean_var_HMSE/(M0-1);
d_mean_var_HMAE = d_mean_var_HMAE/(M0-1);
d_mean_var_QLIKE = d_mean_var_QLIKE/(M0-1);
d_mean_var_R2LOG = d_mean_var_R2LOG/(M0-1);


for i = 1 : M0
    for j = 1 : M0
        temp_Tsq(i,j,1) = (d_mean_MSE(i,j)^2/d_mean_var_MSE(1,i));
        temp_Tsq(i,j,2) = (d_mean_MAE(i,j)^2/d_mean_var_MAE(1,i));
        temp_Tsq(i,j,3) = (d_mean_HMSE(i,j)^2/d_mean_var_HMSE(1,i));
        temp_Tsq(i,j,4) = (d_mean_HMAE(i,j)^2/d_mean_var_HMAE(1,i));
        temp_Tsq(i,j,5) = (d_mean_QLIKE(i,j)^2/d_mean_var_QLIKE(1,i));
        temp_Tsq(i,j,6) = (d_mean_R2LOG(i,j)^2/d_mean_var_R2LOG(1,i));
    end
end

for i = 1 : M0
    Tsq(i,1) = max(temp_Tsq(i,:,1));
    Tsq(i,2) = max(temp_Tsq(i,:,2));
    Tsq(i,3) = max(temp_Tsq(i,:,3));
    Tsq(i,4) = max(temp_Tsq(i,:,4));
    Tsq(i,5) = max(temp_Tsq(i,:,5));
    Tsq(i,6) = max(temp_Tsq(i,:,6));
end