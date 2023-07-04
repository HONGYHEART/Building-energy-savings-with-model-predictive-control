function [A,y] = rLatticeValue(lattice,z)
%在z_linear最后一行添加1
% z_linear = [z_linear;ones(size(z_linear,2),1)'];

z = [z;1];
y = zeros(size(lattice,1),1);
A = zeros(size(lattice,1),size(z,1));

for k = 1:size(lattice,1)
    maxValue = -inf;
    index = [1,1];
    %     k;
    for i = 1:size(lattice{k},2)
        %        i
        minValue = lattice{k}{i}{1} * z;
        minIndex = [i,1];
        for j = 1:size(lattice{k}{i},2)
            if (lattice{k}{i}{j} * z < minValue)
                minValue = lattice{k}{i}{j} * z;
                minIndex = [i,j];
            end
        end
        if(minValue > maxValue)
            maxValue = minValue;
            index = minIndex;
        end
    end
    %    y = [y;maxValue];
    y(k) = maxValue;
    A(k,:) = lattice{k}{index(1)}{index(2)};
%      A(k,:) = index;
end
end