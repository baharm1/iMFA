function [slope] = input_metab_slope(ip_matrix,time_points,missing)

% Fill missing values. Only fills in one consecutive missing value.
if(missing)
    for i = 1:size(ip_matrix,1)
        for j = size(ip_matrix,2)
            if(isnan(ip_matrix(i,j)))
                a1 = ip_matrix(i,j-1);
                a3 = ip_matrix(i,j+1);
                t1 = time_points(j-1);
                t2 = time_points(j);
                t3 = time_points(j+1);
                ip_matrix(i,j) = a1 + (a3-a1)*(t2-t1)/(t3-t1);
            end
        end
    end
end

% Define the slope matrix
slope = zeros(size(ip_matrix,1),numel(time_points)-1);
for j = 1:size(slope,2)
    slope(:,j) = (ip_matrix(:,j+1) - ip_matrix(:,j))/(time_points(j+1)-time_points(j));
end
