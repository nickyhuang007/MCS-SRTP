% function L = h(a,b,k)
clc
clear

Idx = xlsread('C:\Users\Lenovo\Desktop\国创\道琼斯工业平均指数.xls','E2 : E22754');%第一位为0，因为为第一天手动置为0
Idx = Idx';

syms h0 v temp a b k;

temp = h0;
e = sqrt(h0) * v;
for i = 1 : 2
    h1 = (k + b * h0 + a * e^2);
    h0 = h1;
    e = sqrt(h0) * v;
end

L = (-0.5 * log(2*pi) - 0.5 * log(h1) -0.5 * (e^2 / h1));

h0 = temp;

L = subs(L,{h0,v},{1,0.5});
% 
% cal = @(a,b,k)(L);
% 
% pha = mle(Idx,'logpdf',cal,'start',[0,0,0]);

% pha = mle(Idx,'nloglf',L);
% 
% dL_k = diff(L,k);
% dL_b = diff(L,b);
% dL_a = diff(L,a);
% 
% hhh = solve('dL_k = 0','dL_b = 0','dL_a = 0','k','b','a');