%metropolis: many-body interaction
%H=\sum -J2*SiSj - J4*SiSjSmSn - J6*SiSjSmSnSpSq  (近邻+四体+六体)
%T=(2/log(1+sqrt(2)));
%% 清理内存
clc
clear
close all
%% 参数设置  %% 不同温度模拟
n=6;
T=(0:0.1:5);
J2=1;  J4=1;  J6=0;
Thermal=5000;          %弛豫步数
bins=1000;              %bins数目
bsteps=2^8;            %每个bin的步数
% for i=1:length(T)
%     [e(i),m(i),e_error(i),m_error(i),cv(i),ms(i),~,~]=mcmc(n,T(i),Thermal,bins,bsteps,J2,J4,J6);
% end
% figure(1);errorbar(T,e,e_error,'bo:');xlabel('T');ylabel('E');hold on
% figure(2);errorbar(T,m,m_error,'bo:');xlabel('T');ylabel('M');hold on
% figure(3);plot(T,cv,'ro:');xlabel('T');ylabel('C_{v}');hold on
% figure(4);plot(T,ms,'ro:');xlabel('T');ylabel('\chi'); hold on
%% 标度不变（用ms磁化率验证）
% yl=ms*n^(-7/4);
% Tc=(2/log(1+sqrt(2)));
% xl=((T-Tc)/Tc)*n^(1);
% figure(5);plot(xl,yl,'ro:');  xlabel('xl'); ylabel('yl');hold on
%% 计算关联时间
n=[4,8,16,32];
T=(2/log(1+sqrt(2)));
J2=1;  J4=0;  J6=0;
Thermal=5000;          %弛豫步数
bins=1000;              %bins数目
bsteps=1000;            %每个bin的步数
for i=1:length(n)
    [~,~,~,~,~,~,esteps,msteps]=mcmc(n(i),T,Thermal,bins,bsteps,J2,J4,J6);
    [tau(i)]=aaaa(msteps,bins,bsteps);
end
figure;plot(n(1:3),tau,'ro:');  xlabel('L'); ylabel('\tau'); hold on
%%
function [tau]=aaaa(x,bins,bsteps)
j=1;
data_bin=zeros(1,bins);
for i=1:length(x)
    data_bin(j)=data_bin(j)+x(i);
    while i==bsteps*j
        data_bin(j)=data_bin(j)/bsteps;
        j=j+1;
    end
end
tau=1/2 * (bsteps*var(data_bin)/var(x)-1);
end
%% 主函数: 计算当前温度下的物理量
function [e,m,e_error,m_error,cv,ms,e_steps,m_steps]=mcmc(n,T,Thermal,bins,bsteps,J2,J4,J6)
fprintf('temperature is %f\t',T); 
fprintf('L is %f\t',n); 
fprintf('--已开始--\n');
beta=1/T;
[ising]=generate(n);%生成
[nbor]=neighbor(n); %邻居结构
for j=1:Thermal
    [ising]=flip_spin(ising,nbor,n,beta,J2,J4,J6);%一个mc步
end
for j=1:bins*bsteps
    [ising]=flip_spin(ising,nbor,n,beta,J2,J4,J6);%一个mc步
    [e_steps(j),m_steps(j)]=calculate(ising,nbor,n,J2,J4,J6);
end
%求平均
e=mean(e_steps); e_error=error_bar(e_steps,bins,bsteps);
m=mean(m_steps); m_error=error_bar(m_steps,bins,bsteps);
[cv,ms]=calculate_ShMs(e_steps,m_steps,bins,bsteps,T,n);
fprintf('temperature is %f\t',T); fprintf('--已完成--\n');
end
%% 自关联函数
function [avt,ert]=auto_correlation_function(m,ntime,bins,bsteps)
tobs=zeros(ntime,1);
adata=zeros((ntime+1),bins);
i=1; % 时间
for ib =1: bins
    count1 = 1;%判断tobs是否足够包含ntime个时刻的次数
    count2 = 0;%记录平移次数
    aobs = 0;  %平均值
    acor = zeros(ntime,1);%关联
    for ist = 1: bsteps
        if (count1<=ntime)
            tobs (count1)= m(i);
            aobs = aobs + m(i);
        else
            aobs = aobs + m(i);
            count2 = count2 + 1;
            %平移tobs以包含下一时刻的数据
            for z=1:(ntime-1)
                tobs(z) = tobs(z+1);
            end
            tobs(ntime)=m(i);
            %计算不同时刻的关联
            for z = 1: ntime
                acor(z) = acor(z)+ tobs(1)*tobs(z); %计算关联
            end
        end
        count1 = count1 + 1;
        i=i+1;
    end
    aobs=aobs/bsteps;
    acor=acor/count2;
    adata(:,ib)=[aobs',acor']';
end
% 第二步
%Average autocorrelation function
avt=zeros(ntime,1);
for i=1:ntime
    for j=1:bins
        avt(i)=avt(i) + ( adata(i+1,j) - adata(1,j)^2);
    end
end
avt=avt/bins;%对所有bin平均
avt=avt/avt(1);%除以0时刻的值
% Error bars of autocorrelation function
ert=zeros(ntime,1);
abin=zeros(ntime,1);
for j=1:bins
    for i=1:ntime
        abin(i)=adata(i+1,j)-adata(1,j)^2;
    end
    abin=abin/abin(1);
    for i=1:ntime
        ert(i)=ert(i)+(avt(i)-abin(i))^2;
    end
end
ert=sqrt(ert/(bins*(bins-1)));
end
%% 误差
function y=error_bar(x,bins,bsteps)
j=1;
data_bin=zeros(1,bins);
for i=1:length(x)
    data_bin(j)=data_bin(j)+x(i);
    while i==bsteps*j
        data_bin(j)=data_bin(j)/bsteps;
        j=j+1;
    end
end
y=0;
for i=1:bins
    y=y+(data_bin(i)-mean(data_bin))^2;
end
y=sqrt(y/(bins*(bins-1)));
end
%% 初始构型
function [ising]=generate(n)
ising=ones(n*n,1);
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
function [ising]=flip_spin(ising,nbor,n,beta,J2,J4,J6)
for i=1:n*n
    cen=unidrnd(n*n);
    %最近邻
    right=nbor(cen,1);  up=nbor(cen,2);
    left=nbor(cen,3);   down=nbor(cen,4);
    
    %次近邻
    rightup=nbor(right,2);       leftup=nbor(left,2);
    rightdown=nbor(right,4);     leftdown=nbor(left,4);
    
    %六体
    upup=nbor(up,2);downdown=nbor(down,4);upupleft=nbor(upup,3);
    upupright=nbor(upup,1);downdownleft=nbor(downdown,3);downdownright=nbor(downdown,1);
    
    deE=2 * J2 * ising(cen) * ( ising(right) + ising(left) + ising(up) + ising(down) ) ;
    
    deE=deE + 2 * J4 * ising(cen) * ( ising(up)*ising(rightup)*ising(right) +ising(up)*ising(leftup)*ising(left) +ising(down)*ising(rightdown)*ising(right) +ising(down)*ising(leftdown)*ising(left) ) ;
    
    s1=ising(up)*ising(upup)*ising(upupright)*ising(rightup)*ising(right) + ising(up)*ising(rightup)*ising(right)*ising(rightdown)*ising(down);
    s2=ising(down)*ising(downdown)*ising(downdownright)*ising(rightdown)*ising(right) + ising(up)*ising(upup)*ising(upupleft)*ising(leftup)*ising(left);
    s3=ising(up)*ising(leftup)*ising(left)*ising(leftdown)*ising(down) + ising(down)*ising(downdown)*ising(downdownleft)*ising(leftdown)*ising(left);
    deE=deE + 2 * J6 * ising(cen) * ( s1 + s2 + s3 );
    
    if rand<exp(-deE*beta)
        ising(cen)=-ising(cen);
    end
end
end
%% 计算参量
function [e,m]=calculate(ising,nbor,n,J2,J4,J6)
e=0.0;
for i=1:n*n
    e=e - J2 * ising(i) * ( ising(nbor(i,1)) + ising(nbor(i,2)) );
    e=e - J4 * ising(i) * ising(nbor(i,1)) * ising(nbor(i,2)) * ising(nbor(nbor(i,2),1));
    e=e - J6 * ( ising(i) * ising(nbor(i,1)) * ising(nbor(i,2)) * ising(nbor(nbor(i,1),2)) * ising(nbor(nbor(i,2),2)) * ising(nbor(nbor(nbor(i,2),2),1)) );
end
e=e/n^2;
m=abs(sum(ising))/n^2;
end
%% 计算比热、磁化率
function [cv,ms]=calculate_ShMs(e,m,bins,bsteps,T,n)
%---------------不用block平均-------------------
%cv=(mean(e.^2)-mean(e)^2)/T^2;
%ms=(mean(m.^2)-mean(m)^2)/T;
%---------------用block平均---------------------
j=1;
e_bins=zeros(1,bins);
m_bins=zeros(1,bins);
e2_bins=zeros(1,bins);
m2_bins=zeros(1,bins);
for i=1:length(e)
    e_bins(j)=e_bins(j)+e(i);
    m_bins(j)=m_bins(j)+m(i);
    e2_bins(j)=e2_bins(j)+e(i)^2;
    m2_bins(j)=m2_bins(j)+m(i)^2;
    while i==bsteps*j
        e_bins(j)=e_bins(j)/bsteps;
        m_bins(j)=m_bins(j)/bsteps;
        e2_bins(j)=e2_bins(j)/bsteps;
        m2_bins(j)=m2_bins(j)/bsteps;
        j=j+1;
    end
end
cv=n^2*(mean(e2_bins)-mean(e_bins)^2)/T^2;
ms=n^2*(mean(m2_bins)-mean(m_bins)^2)/T;
end