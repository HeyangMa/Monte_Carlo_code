%穷举  H=\sum -J2*SiSj
% clear;clc;
% close all;
%% set
n=3;
J2=1;
T_list= (1:0.1:5);
h_list= 0;%(-0.1:0.01:0.1)
pbc=1;
fbc=0;
%% prepare
nbor = neighbor_pbc(n);
e = [];
m = [];
Cv = [];
ms = [];
free_energy = [];
%% main
for j=1:length(T_list)
    for i=1:length(h_list)
        beta=1/T_list(j);
        h=h_list(i);
        [e(end+1),m(end+1),Cv(end+1),ms(end+1),free_energy(end+1)]=main_func(n,pbc,fbc,nbor,beta,J2,h);
        fprintf('beta = %f  ', beta);
        fprintf('h = %f\n',   h);
    end
end
%% figure
if length(T_list)~=1; x=T_list; else; x=h_list;end
figure(1);hold on;plot(x, e,'k*-');ylabel('Energy')
figure(2);hold on;plot(x, m,'k*-');ylabel('mag')
figure(3);hold on;plot(x,Cv,'k*-');ylabel('Cv')
figure(4);hold on;plot(x,ms,'k*-');ylabel('ms')
figure(5);hold on;plot(x,free_energy,'k*-');ylabel('free energy')

%% main function
function [energy_per_site,mag_per_site,Cv,ms,free_energy_per_site]=main_func(n,pbc,fbc,nbor,beta,J2,h)
% 数组大小
z=[];
ising = -ones(1,n*n);
energy = zeros(2^(n*n),1);
m = zeros(2^(n*n),1);

% 初始全为-1,计算
[energy(1)] = calculate(ising, nbor, n, J2, h, pbc, fbc);
m(1)=abs(sum(ising));
z(1) = exp(-beta*energy(1));
% 2^n*n-1个值
for i = 1:length(energy)-1
    a = str2double(dec2bin(i));
    a = num2str(a, append('%0', num2str(n^2), 'd'));
    for j=1:n*n
        if str2double(a(j)) == 1
            ising(j) = 1;
        else
            ising(j) = -1;
        end
    end
    [energy(i+1)] = calculate(ising, nbor, n, J2, h, pbc, fbc);
    m(i+1)=abs(sum(ising));
    z(i+1) = exp(-beta*energy(i+1));
end
a=unique(energy);
Z = sum(z);                                                       % 求配分函数
Eng = sum(energy .* exp(-beta * energy)) / Z;                     % 计算平均能量
Eng2 = sum(energy .* energy .* exp(-beta * energy)) / Z;          % 计算平均能量的平方值
mag = sum(m .* exp(-beta * energy)) / Z;                          % 计算平均磁化强度
mag2 = sum(m .* m .* exp(-beta * energy)) / Z;                    % 计算平均磁化强度的平方值

energy_per_site = Eng/n^2;
mag_per_site = mag/n^2;
Cv = (Eng2 - Eng^2) * beta^2 / n^2;          % 计算比热
ms = (mag2 - mag^2) * beta^1 / n^2;          % 计算磁化率
free_energy_per_site = -log(Z)/beta/n^2;     % 计算自由能
end
%% generate neighbor
function [nbor]=neighbor_pbc(n)
nbor=zeros(n*n,4);
for ispin=1:n*n
    iy=fix((ispin-1)/n)+1;
    ix=ispin-(iy-1)*n;
    ixp=ix+1-fix(ix/n)*n;       %右侧点, x坐标
    iyp=iy+1-fix(iy/n)*n;       %上侧点, y坐标
    ixm=ix-1+fix((n-ix+1)/n)*n; %左侧点, x坐标
    iym=iy-1+fix((n-iy+1)/n)*n; %下侧点, y坐标
    nbor(ispin,1)=(iy-1)*n+ixp;%右邻居
    nbor(ispin,2)=(iyp-1)*n+ix;%上邻居
    nbor(ispin,3)=(iy-1)*n+ixm;%左邻居
    nbor(ispin,4)=(iym-1)*n+ix;%下邻居
end
end
%% calculate_energy
function [energy]=calculate(ising,nbor,n,J2,h,pbc,fbc)
%pbc
if pbc==1
    energy=0.0;
    for i=1:n*n
        energy=energy - J2 * ising(i) * ( ising(nbor(i,1)) + ising(nbor(i,2)) );
    end
    for i=1:n*n
        energy=energy - h * ising(i);
    end
end

%fbc
if fbc==1
    energy=0.0;
    for i=1:n*n
        iy=fix((i-1)/n)+1;
        ix=i-(iy-1)*n;
        if ix~=n && iy~=n
            energy=energy - J2 * ising(i) * ( ising(nbor(i,1)) + ising(nbor(i,2)) );
        end
        if ix==n && iy~=n
            energy=energy - J2 * ising(i) * ( ising(nbor(i,2)) );
        end
        if iy==n && ix~=n
            energy=energy - J2 * ising(i) * ( ising(nbor(i,1)) );
        end
    end
    for i=1:n*n
        energy=energy - h * ising(i);
    end
end
end