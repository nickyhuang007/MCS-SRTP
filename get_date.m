function [y,m,d] = get_date(date)
p=strfind(date,'/');
N = size(date,1);
y = zeros(N,1);
m = zeros(N,1);
d = zeros(N,1); 
for i=1:N
    y(i,1)=str2double(date{i}(1:p{i}(1)-1));
    m(i,1)=str2double(date{i}(p{i}(1)+1:p{i}(2)-1));
    d(i,1)=str2double(date{i}(p{i}(2)+1:end));
end