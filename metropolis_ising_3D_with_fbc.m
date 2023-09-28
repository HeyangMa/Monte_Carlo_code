%metropolis: 3D ising Nearest-Neighbor interaction with fbc
%H=\sum -J2*SiSj
%% 清理内存
clc
clear
close all
%% 参数设置
n=3;
T=(0:0.1:4);
J2=1;
Thermal=1000;           %弛豫步数
bins=500;              %bins数目
bsteps=10;            %每个bin的步数
%% 主函数
[ising_initial]=generate(n);%生成
[nbor]=neighbor(n); %邻居结构
for i=1:length(T)
    beta=1/T(i);
    ising = ising_initial;
    for j=1:Thermal
        [ising]=flip_spin(ising,nbor,n,beta,J2);%一个mc步
    end
    for j=1:bins*bsteps
        [ising]=flip_spin(ising,nbor,n,beta,J2);%一个mc步
        [e_steps(j),m_steps(j)]=calculate(ising,nbor,n,J2);
    end
    %求平均
    e(i)=mean(e_steps);
    m(i)=mean(m_steps);
    %cv, ms
    cv(i) = beta^2 * var( n^3 * e_steps);
    ms(i) = beta * n^3 * var(m_steps);
    %binder
    R2 = mean(m_steps.^4) / mean(m_steps.^2)^2;
    B_r(i) = (3 / 2) * (1 - R2);
    fprintf('temperature is %f\t',T(i)); fprintf('--已完成--\n');
end
%% plot
figure(1);plot(T,e,'rs-');  xlabel('T');ylabel('E');    hold on
figure(2);plot(T,m,'rs-');  xlabel('T');ylabel('M');    hold on
figure(3);plot(T,cv,'rs-'); xlabel('T');ylabel('C_{v}');hold on
figure(4);plot(T,ms,'rs-'); xlabel('T');ylabel('\chi'); hold on
figure(5);plot(T,B_r,'rs-');xlabel('T');ylabel('B_{r}');hold on

%% 初始构型
function [ising]=generate(n)
ising=ones(n^3,1);
end
%% 邻居表
function [neigh]=neighbor(n)
%generate table of near neighbor
neigh = zeros( n^3, 6);
for ispin = 1:n^3
    iz = fix((ispin-1)/n^2) + 1;
    iy = fix((ispin-((iz-1)*n^2)-1)/n) + 1;
    ix = (ispin-(iz-1)*n^2) - (iy-1) * n;
    
    ixp = ix + 1 - fix(ix/n) * n;
    iyp = iy + 1 - fix(iy/n) * n;
    izp = iz + 1 - fix(iz/n) * n;
    ixm = ix - 1 + fix((n-ix+1)/n) * n;
    iym = iy - 1 + fix((n-iy+1)/n) * n;
    izm = iz - 1 + fix((n-iz+1)/n) * n;
    neigh(ispin,1) = (iz-1)  * n^2 + (iy-1)  * n + ixp; %右邻居
    neigh(ispin,2) = (iz-1)  * n^2 + (iyp-1) * n + ix;  %上邻居
    neigh(ispin,3) = (iz-1)  * n^2 + (iy-1)  * n + ixm; %左邻居
    neigh(ispin,4) = (iz-1)  * n^2 + (iym-1) * n + ix;  %下邻居
    neigh(ispin,5) = (izp-1) * n^2 + (iy-1)  * n + ix;  %上层邻居
    neigh(ispin,6) = (izm-1) * n^2 + (iy-1)  * n + ix;  %下层邻居
    %fbc
%     if ix==1
%         neigh(ispin,3)=0;
%     elseif ix==n
%         neigh(ispin,1)=0;
%     end
%     if iy==1
%         neigh(ispin,4)=0;
%     elseif iy==n
%         neigh(ispin,2)=0;
%     end
%     if iz==1
%         neigh(ispin,6)=0;
%     elseif iz==n
%         neigh(ispin,5)=0;
%     end
end
end
%% 一步mc
function [ising]=flip_spin(ising,nbor,n,beta,J2)
for i=1:n^3
    cen=unidrnd(n^3);
    %最近邻
    if nbor(cen,1)~=0
        right  = ising(nbor(cen,1));
    else
        right  = 0;
    end
    if nbor(cen,2)~=0
        up     = ising(nbor(cen,2));
    else
        up     = 0;
    end
    if nbor(cen,3)~=0
        left   = ising(nbor(cen,3));
    else
        left   = 0;
    end
    if nbor(cen,4)~=0
        down   = ising(nbor(cen,4));
    else
        down   = 0;
    end
    if nbor(cen,5)~=0
        top    = ising(nbor(cen,5));
    else
        top    = 0;
    end
    if nbor(cen,6)~=0
        bottom = ising(nbor(cen,6));
    else
        bottom = 0;
    end
    deE=2 * J2 * ising(cen) * ( right + left + up + down + top + bottom ) ;
    if rand<exp(-deE*beta)
        ising(cen)=-ising(cen);
    end
end
end
%% 计算参量
function [e,m]=calculate(ising,nbor,n,J2)
e=0.0;
for i=1:n^3
    if nbor(i,1)~=0
        right  = ising(nbor(i,1));
    else
        right  = 0;
    end
    if nbor(i,2)~=0
        up     = ising(nbor(i,2));
    else
        up     = 0;
    end
    if nbor(i,5)~=0
        top    = ising(nbor(i,5));
    else
        top    = 0;
    end
    e = e - J2 * ising(i) * ( right + up + top );
end
e=e/n^3;
m=abs(sum(ising))/n^3;
end