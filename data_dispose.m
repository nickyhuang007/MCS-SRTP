clc             
clear

Idx = xlsread('C:\Users\Lenovo\Desktop\国创\000001.SS.csv','B2 : E7253');%第一位为0，因为为第一天手动置为0
[~,Trd] = xlsread('C:\Users\Lenovo\Desktop\国创\000001.SS.csv','A2 : A7253');


count = 0;
for i = 1 : length(Idx)
    if(isnan(Idx(i,3)))
        count = count + 1;
    end
end

pos = zeros(1,count);
new_Idx = zeros(length(Idx)-count,4);

for k = 1 : 4
    j = 1;
    h = 1;
  for i = 1 : length(Idx)
      if(isnan(Idx(i,k)))
          pos(h) = i;
          h = h + 1;
      else
          new_Idx(j,k) = Idx(i,k);
          j = j + 1;
      end
  end
end

new_Trd = cell(length(Idx) - count,1);
j = 1;h = 1;
for i = 1 : length(Idx)
    if(i == pos(1,h))
        h = h + 1;
        if(h > count)
            h = count;
        end
    else
        new_Trd{j,1} = Trd{i,1};
        j = j + 1;
    end
end

xlswrite('C:\Users\Lenovo\Desktop\国创\new_000001.SS.xlsx',new_Idx,'B2 : E7238');
xlswrite('C:\Users\Lenovo\Desktop\国创\new_000001.SS.xlsx',new_Trd,'A2 : A7238');


%[year,month,day] = get_date(Trd);

