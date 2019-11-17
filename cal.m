clc
clear

%Idx = xlsread('C:\Users\Lenovo\Desktop\国创\道琼斯工业平均指数.xls','E3 : E22754');%第一位为0，因为为第一天手动置为0
 
Idx = xlsread('C:\Users\Lenovo\Desktop\国创\new_000001.SS.xlsx','E2 : E7238');%多取一个数据来求第一位的方差

H = 7000; %估计样本数
M = length(Idx) - 7000; %测试样本数
M0 = 3;  %要估计的模型数
est_Idx = Idx(1:H);
test_Idx = Idx(H+1:length(Idx));

com_data = zeros(M,M0);%com值估计后的值
com_arch = zeros(M,M0);
com_mean = zeros(M,M0);
com_diff = zeros(M,M0);
com_SSE = zeros(M,M0);


[com_data(:,1),com_arch(:,1),com_mean(:,1),com_diff(:,1),com_SSE(:,1)] = cal_garch(test_Idx,1.52284e-05,0.765176,0.234824,-0.001755);
[com_data(:,2),com_arch(:,2),com_mean(:,2),com_diff(:,2),com_SSE(:,2)] = cal_gjr_garch(test_Idx,1.01085e-05,0.812004,0.227273,-0.0785537,-0.001755);
[com_data(:,3),com_arch(:,3),com_mean(:,3),com_diff(:,3),com_SSE(:,3)] = cal_egarch(test_Idx,-0.0561225,0.989209,0.22702,0.0384927,-0.001755);
% a_diff = 1/a_diff;
% b_diff = 1/b_diff;
% c_diff = 1/c_diff;
% 
% A = [a_mean,a_diff,a_SSE;
%      b_mean,b_diff,b_SSE;
%      c_mean,c_diff,c_SSE];
% W = [1/3,1/3,1/3];
% hhh = TOPSIS(A,W);
% Idx = xlsread('C:\Users\Lenovo\Desktop\国创\new_000001.SS.xlsx','E7001 : E7238');
data_arch = zeros(length(Idx),1);

for i = 2 : length(Idx)
    data_arch(i,1) = ((Idx(i,1) - Idx(i-1,1))/Idx(i-1,1))^2;    
end

MSE = zeros(M,M0);
MAE = zeros(M,M0);
HMSE = zeros(M,M0);
HMAE = zeros(M,M0);
QLIKE = zeros(M,M0);
R2LOG = zeros(M,M0);

for k = 1 : M0
  for j = 1 : M
    for i = 1 : j
        MSE(j,k) = MSE(j,k) + (data_arch(H+i) - com_arch(i,k))^2;
        MAE(j,k) = MAE(j,k) + abs(data_arch(H+i) - com_arch(i,k));
        HMSE(j,k) = HMSE(j,k) + (1 - (com_arch(i,k)/data_arch(H+i)))^2;
        HMAE(j,k) = HMAE(j,k) + abs(1 - (com_arch(i,k)/data_arch(H+i)));
        QLIKE(j,k) = QLIKE(j,k) + log(com_arch(i,k)) + data_arch(H+i)/com_arch(i,k);
        R2LOG(j,k) = R2LOG(j,k) + (log(data_arch(H+i))/com_arch(i,k))^2;
    end
  end
end

for j = 1 : M0
  for i = 1 : M
    MSE(i,j) = MSE(i,j)/i;
    MAE(i,j) = MAE(i,j)/i;
    HMSE(i,j) = HMSE(i,j)/i;
    HMAE(i,j) = HMAE(i,j)/i;
    QLIKE(i,j) = QLIKE(i,j)/i;
    R2LOG(i,j) = R2LOG(i,j)/i;
  end
end

d_MSE = zeros(M0,M0,length(MSE)); %记录3个模型分别于另外两个模型之间的差,分6种函数
d_MAE = zeros(M0,M0,length(MAE));
d_HMSE = zeros(M0,M0,length(HMSE));
d_HMAE = zeros(M0,M0,length(HMAE));
d_QLIKE = zeros(M0,M0,length(QLIKE));
d_R2LOG = zeros(M0,M0,length(R2LOG));

for i = 1 : M0
    for j = 1 : M0
        if(i ~= j)
            d_MSE(i,j,:) = MSE(:,i) - MSE(:,j);
            d_MAE(i,j,:) = MAE(:,i) - MAE(:,j);
            d_HMSE(i,j,:) = HMSE(:,i) - HMSE(:,j);
            d_HMAE(i,j,:) = HMAE(:,i) - HMAE(:,j);
            d_QLIKE(i,j,:) = QLIKE(:,i) - QLIKE(:,j);
            d_R2LOG(i,j,:) = R2LOG(:,i) - R2LOG(:,j);
        end
    end
end

ran_d_MSE = zeros(M0,M0,length(MSE)); %记录3个模型分别于另外两个模型之间的差,分6种函数
ran_d_MAE = zeros(M0,M0,length(MAE));
ran_d_HMSE = zeros(M0,M0,length(HMSE));
ran_d_HMAE = zeros(M0,M0,length(HMAE));
ran_d_QLIKE = zeros(M0,M0,length(QLIKE));
ran_d_R2LOG = zeros(M0,M0,length(R2LOG));

N = 1000;  %Brootstrap迭代次数
Tr = zeros(N,M0,6);
Tsq = zeros(N,M0,6);

for n = 1 : N
  for i = 1 : 3
      for j = 1 : 3
          ran = randi([1 M],1,M);
          for k = 1 : M
              ran_d_MSE(i,j,k) = d_MSE(i,j,ran(k));
              ran_d_MAE(i,j,k) = d_MAE(i,j,ran(k));
              ran_d_HMSE(i,j,k) = d_HMSE(i,j,ran(k));
              ran_d_HMAE(i,j,k) = d_HMAE(i,j,ran(k));
              ran_d_QLIKE(i,j,k) = d_QLIKE(i,j,ran(k));
              ran_d_R2LOG(i,j,k) = d_R2LOG(i,j,ran(k));
          end
      end
  end
  [Tr(n,:,:),Tsq(n,:,:)] = get_statistics(ran_d_MSE,ran_d_MAE,ran_d_HMSE,ran_d_HMAE,ran_d_QLIKE,ran_d_R2LOG,M0,M);
end

for i = 1 : M0
    for j = 1 : 6
        Tr(:,i,j) = sort(Tr(:,i,j));
        Tsq(:,i,j) = sort(Tsq(:,i,j));
    end
end

  [real_Tr,real_Tsq] = get_statistics(d_MSE,d_MAE,d_HMSE,d_HMAE,d_QLIKE,d_R2LOG,M0,M);
  
  pos_real_Tr = zeros(M0,6);
  pos_real_Tsq = zeros(M0,6);
  
for i = 1 : M0
    for j = 1 : 6
        for k = 1 : N-1
            if((real_Tr(i,j) >= Tr(k,i,j)) && (real_Tr(i,j) <= Tr(k+1,i,j)))
                pos_real_Tr(i,j) = k + (real_Tr(i,j) - Tr(k,i,j))/(Tr(k+1,i,j) - Tr(k,i,j));
            end
            
            if((real_Tsq(i,j) >= Tsq(k,i,j)) && (real_Tsq(i,j) <= Tsq(k+1,i,j)))
                pos_real_Tsq(i,j) = k + (real_Tsq(i,j) - Tsq(k,i,j))/(Tsq(k+1,i,j) - Tsq(k,i,j));
            end
        end
    end
end
  
p_Tr = zeros(M0,6);
p_Tsq = zeros(M0,6);

for i = 1 : M0
    for j = 1 : 6
        p_Tr(i,j) = 1 - abs(1 - (2 * pos_real_Tr(i,j))/N);
        p_Tsq(i,j) = 1 - abs(1 - (2 * pos_real_Tsq(i,j))/N);
    end
end


% d_mean_MSE = zeros(M0,M0);
% d_mean_MAE = zeros(M0,M0);
% d_mean_HMSE = zeros(M0,M0);
% d_mean_HMAE = zeros(M0,M0);
% d_mean_QLIKE = zeros(M0,M0);
% d_mean_R2LOG = zeros(M0,M0);
% 
% for i = 1 : M0
%     for j = 1 : M0
%         for k = 1 : M
%             d_mean_MSE(i,j) = d_mean_MSE(i,j) + d_MSE(i,j,k);
%             d_mean_MAE(i,j) = d_mean_MAE(i,j) + d_MAE(i,j,k);
%             d_mean_HMSE(i,j) = d_mean_HMSE(i,j) + d_HMSE(i,j,k);
%             d_mean_HMAE(i,j) = d_mean_HMAE(i,j) + d_HMAE(i,j,k);
%             d_mean_QLIKE(i,j) = d_mean_QLIKE(i,j) + d_QLIKE(i,j,k);
%             d_mean_R2LOG(i,j) = d_mean_R2LOG(i,j) + d_R2LOG(i,j,k);
%         end
%     end
% end
% 
% d_mean_MSE = d_mean_MSE/M;
% d_mean_MAE = d_mean_MAE/M;
% d_mean_HMSE = d_mean_HMSE/M;
% d_mean_HMAE = d_mean_HMAE/M;
% d_mean_QLIKE = d_mean_QLIKE/M;
% d_mean_R2LOG = d_mean_R2LOG/M;
% 
% d_var_MSE = zeros(M0,M0);
% d_var_MAE = zeros(M0,M0);
% d_var_HMSE = zeros(M0,M0);
% d_var_HMAE = zeros(M0,M0);
% d_var_QLIKE = zeros(M0,M0);
% d_var_R2LOG = zeros(M0,M0);
% 
% for i = 1 : M0
%     for j = 1 : M0
%         for k = 1 : M
%             d_var_MSE(i,j) = d_var_MSE(i,j) + (d_MSE(i,j,k) - d_mean_MSE(i,j))^2;
%             d_var_MAE(i,j) = d_var_MAE(i,j) + (d_MAE(i,j,k) - d_mean_MAE(i,j))^2;
%             d_var_HMSE(i,j) = d_var_HMSE(i,j) + (d_HMSE(i,j,k) - d_mean_HMSE(i,j))^2;
%             d_var_HMAE(i,j) = d_var_HMAE(i,j) + (d_HMAE(i,j,k) - d_mean_HMAE(i,j))^2;
%             d_var_QLIKE(i,j) = d_var_QLIKE(i,j) + (d_QLIKE(i,j,k) - d_mean_QLIKE(i,j))^2;
%             d_var_R2LOG(i,j) = d_var_R2LOG(i,j) + (d_R2LOG(i,j,k) - d_mean_R2LOG(i,j))^2;
%         end
%     end
% end
% 
% d_var_MSE = d_var_MSE/M;
% d_var_MAE = d_var_MAE/M;
% d_var_HMSE = d_var_HMSE/M;
% d_var_HMAE = d_var_HMAE/M;
% d_var_QLIKE = d_var_QLIKE/M;
% d_var_R2LOG = d_var_R2LOG/M;
% 
% Tr = zeros(M0,6);
% Tsq = zeros(M0,6);
% 
% temp_Tr = zeros(M0,M0,6);
% temp_Tsq = zeros(M0,M0,6);
% 
% for i = 1 : M0
%     for j = 1 : M0
%         temp_Tr(i,j,1) = (abs(d_mean_MSE(i,j))/sqrt(d_var_MSE(i,j)));
%         temp_Tr(i,j,2) = (abs(d_mean_MAE(i,j))/sqrt(d_var_MAE(i,j)));
%         temp_Tr(i,j,3) = (abs(d_mean_HMSE(i,j))/sqrt(d_var_HMSE(i,j)));
%         temp_Tr(i,j,4) = (abs(d_mean_HMAE(i,j))/sqrt(d_var_HMAE(i,j)));
%         temp_Tr(i,j,5) = (abs(d_mean_QLIKE(i,j))/sqrt(d_var_QLIKE(i,j)));
%         temp_Tr(i,j,6) = (abs(d_mean_R2LOG(i,j))/sqrt(d_var_R2LOG(i,j)));
%     end
% end
% 
% for i = 1 : M0
%     Tr(i,1) = max(temp_Tr(i,:,1));
%     Tr(i,2) = max(temp_Tr(i,:,2));
%     Tr(i,3) = max(temp_Tr(i,:,3));
%     Tr(i,4) = max(temp_Tr(i,:,4));
%     Tr(i,5) = max(temp_Tr(i,:,5));
%     Tr(i,6) = max(temp_Tr(i,:,6));
% end
% 
% for i = 1 : M0
%     for j = 1 : M0
%         temp_Tsq(i,j,1) = (d_mean_MSE(i,j)^2/d_var_MSE(i,j));
%         temp_Tsq(i,j,2) = (d_mean_MAE(i,j)^2/d_var_MAE(i,j));
%         temp_Tsq(i,j,3) = (d_mean_HMSE(i,j)^2/d_var_HMSE(i,j));
%         temp_Tsq(i,j,4) = (d_mean_HMAE(i,j)^2/d_var_HMAE(i,j));
%         temp_Tsq(i,j,5) = (d_mean_QLIKE(i,j)^2/d_var_QLIKE(i,j));
%         temp_Tsq(i,j,6) = (d_mean_R2LOG(i,j)^2/d_var_R2LOG(i,j));
%     end
% end
% 
% for i = 1 : M0
%     Tsq(i,1) = max(temp_Tsq(i,:,1));
%     Tsq(i,2) = max(temp_Tsq(i,:,2));
%     Tsq(i,3) = max(temp_Tsq(i,:,3));
%     Tsq(i,4) = max(temp_Tsq(i,:,4));
%     Tsq(i,5) = max(temp_Tsq(i,:,5));
%     Tsq(i,6) = max(temp_Tsq(i,:,6));
% end
% x = 1 : length(com_arch);
% plot(x,com_data(:,1),x,com_data(:,2),x,com_data(:,3));
% title('收益率的变化曲线');
% xlabel('天数');
% ylabel('收益率');
% legend('Garch','GJR-Garch','e-garch');
%plot(x,new_Idx(1,:),'b',x,data(1,:),'g');