function [mean_d,std_d,data_norm] = normalized_fun(data_matrix)

mean_d = mean(data_matrix);
std_d = std(data_matrix);
dNorm_x = size(data_matrix,2);
data_norm = zeros(size(data_matrix,1),size(data_matrix,2));
for i = 1:1:dNorm_x
    data_norm(:,i) = (data_matrix(:,i)-mean_d(i))/std_d(i);
end

end

