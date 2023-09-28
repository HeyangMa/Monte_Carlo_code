%metropolis: 2D ising Nearest-Neighbor interaction
%H=\sum -J2*SiSj
%T=(2/log(1+sqrt(2)));
%% 清理内存
clc;clear;
close all
%% 参数设置
n=4;
T=(0:0.1:3);%T=(2/log(1+sqrt(2)));
J2=1;
Thermal=10000;          %弛豫步数
bins=10;               %bins数目
bsteps=200;            %每个bin的步数
%% 主函数
for i=1:length(T)
    beta=1/T(i);
    [ising]=generate(n);%生成
    [nbor]=neighbor(n); %邻居结构
    for j=1:Thermal
        [ising]=flip_spin(ising,nbor,n,beta,J2);%一个mc步
%         [ising]=geometric_update(ising,nbor,n,beta,J2);
    end
    for j=1:bins*bsteps
        [ising]=flip_spin(ising,nbor,n,beta,J2);%一个mc步
%         [ising]=geometric_update(ising,nbor,n,beta,J2);
        [e_steps(j),m_steps(j)]=calculate(ising,nbor,n,J2);
    end
    %求平均
    e(i)=mean(e_steps); 
    m(i)=mean(m_steps);
    cv(i)=beta^2 * var(e_steps);
    ms(i)=beta * var(n^2 * m_steps);
    fprintf('temperature is %f\t',T(i)); fprintf('--已完成--\n');
end
figure(1);plot(T,e,'b*:');xlabel('T');ylabel('E');hold on
figure(2);plot(T,m,'b*:');xlabel('T');ylabel('M');hold on
figure(3);plot(T,cv,'b*:');xlabel('T');ylabel('E');hold on
figure(4);plot(T,m,'b*:');xlabel('T');ylabel('M');hold on

%% 初始构型
function [ising]=generate(n)
ising=ones(n*n,1);
% for i=1:n^2
%     if rand>0.5
%         ising(i)=-1;
%     end
% end
end
%% 邻居表
function [nbor]=neighbor(n)
%generate table of near neighbor
nbor=zeros(n*n,4);
for ispin=1:n*n
    iy= fix((ispin-1)/n)+1;
    ix=ispin-(iy-1)*n;
    ixp=ix+1-fix(ix/n)*n;
    iyp=iy+1-fix(iy/n)*n;
    ixm=ix-1+fix((n-ix+1)/n)*n;
    iym=iy-1+fix((n-iy+1)/n)*n;
    nbor(ispin,1)=(iy-1)*n+ixp;%右邻居
    nbor(ispin,2)=(iym-1)*n+ix;%上邻居
    nbor(ispin,3)=(iy-1)*n+ixm;%左邻居
    nbor(ispin,4)=(iyp-1)*n+ix;%下邻居
end
end
%% 一步mc
function [ising]=flip_spin(ising,nbor,n,beta,J2)
for i=1:n*n
    cen=unidrnd(n*n);
    %最近邻
    right=nbor(cen,1);  up=nbor(cen,2);
    left=nbor(cen,3);   down=nbor(cen,4);
    
    deE=2 * J2 * ising(cen) * ( ising(right) + ising(left) + ising(up) + ising(down) ) ;
    
    if rand<exp(-deE*beta)
        ising(cen)=-ising(cen);
    end
end
end
%% 计算参量
function [e,m]=calculate(ising,nbor,n,J2)
e=0.0;
for i=1:n*n
    e=e - J2 * ising(i) * ( ising(nbor(i,1)) + ising(nbor(i,2)) );
end
e=e/n^2;
m=abs(sum(ising))/n^2;
end
%% GCA
function [ising]=geometric_update(ising,nbor,n,beta,J2)
%初始点
cen = unidrnd(n*n);
%反射点
ref = unidrnd(n*n);
%exchange
ispin = ising(cen);
ising(cen) = ising(ref);
ising(ref) = ispin;
%堆栈
cluster = zeros(n^2, 1);
flag = zeros(n^2, 1);
rlag = zeros(n^2, 1);
%表示把current_site放入了需要查找的队列
flag(1) = cen;
rlag(1) = ref;
%修改值为1，表示该格点已经被选择过
cluster(cen) = 1;
cluster(ref) = 1;
%表示待查找的点的个数+1，如果为1就结束
counter = 2;
%开始查找
while counter > 1
    %从最后一个被加入的点开始查找
    counter = counter - 1;
    cen = flag(counter);
    ref = rlag(counter);
    flag(counter) = 0;
    rlag(counter) = 0;
    %-----%
    for i = 1:4 %J2
        %查找每一个近邻
        next_site = nbor(cen, i);
        if i == 1
            next_site_ref = nbor(ref, 3);
        elseif i == 2
            next_site_ref = nbor(ref, 4);
        elseif i == 3
            next_site_ref = nbor(ref, 1);
        elseif i == 4
            next_site_ref = nbor(ref, 2);
        end
        %需要满足条件
        E_old = ising(cen)*ising(next_site) + ising(ref)*ising(next_site_ref);
        E_new = ising(cen)*ising(next_site_ref) + ising(ref)*ising(next_site);
        deltaE = -J2 * ( E_old - E_new );
        %P=e^(- beta * delta{E}) therefore,
        probability = 1 - exp(-1 * beta * deltaE);
        if deltaE > 0 && rand<probability && cluster(next_site) == 0
            ispin = ising(next_site);
            ising(next_site) = ising(next_site_ref);
            ising(next_site_ref) = ispin;
            %put two sites into one cluster
            cluster(next_site) = 1;
            cluster(next_site_ref) = 1;
            %将满足条件的点放入待查找的队列中
            flag(counter) = next_site;
            rlag(counter) = next_site_ref;
            counter = counter + 1;
        end
    end
    %----------------------
end
end