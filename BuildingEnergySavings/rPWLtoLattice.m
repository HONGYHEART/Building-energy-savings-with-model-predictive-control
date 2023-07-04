function [lattice_f,S] = rPWLtoLattice(PWLs,z_linear)
z_linear = [z_linear;ones(1,size(z_linear,2))];
%PWLsΪ������Գ�ƽ�湹��,z_linearΪ��Ӧ�����Ի���
N_linear_function = size(PWLs{1},1);
%���崢��lattice��Ԫ��,�ṹ����S
lattices = {};
S = {};
for k = 1:N_linear_function
    k;
    %����ÿ������ƽ���lattice,�ṹ������ʱ����2 S_temp2
    lattice = {};
    S_temp2 = [];
    for i = 1:size(PWLs,2);%��ʾ���Ի���
        lattice{i}{1} = PWLs{i}(k,:);
        S_temp1 = [];
        for j = 1:size(PWLs,2)
            if(i==j)
                lattice{i} = {lattice{i}{:},PWLs{j}(k,:)};
                S_temp1 = [S_temp1,1];
            elseif(PWLs{j}(k,:)*z_linear(:,i) > PWLs{i}(k,:)*z_linear(:,i))
                lattice{i} = {lattice{i}{:},PWLs{j}(k,:)};
                S_temp1 = [S_temp1,1];
            else
                S_temp1 = [S_temp1,0];
            end
        end
        S_temp2 = [S_temp2;S_temp1];
        lattice{i} = {lattice{i}{2:end}};
    end
    lattices = [lattices;lattice];
    S = [S;{S_temp2}];
end
for i = 1:N_linear_function
    lattice_f{i,1} = {lattices{i,:}};
end
end