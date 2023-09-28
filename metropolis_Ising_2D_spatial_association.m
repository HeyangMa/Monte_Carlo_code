%metropolis: many-body interaction
%H=\sum -J2*SiSj - J4*SiSjSmSn - J6*SiSjSmSnSpSq  (近邻+四体+六体)
%T=(2/log(1+sqrt(2)));
%% 清理内存
clc
clear
close all
%% 参数设置
n=10;
T=4;
J2=1;
Thermal=1000;          %弛豫步数
bins=100;              %bins数目
bsteps=100;            %每个bin的步数
%% 不同温度模拟
[S_array]=mcmc(n,T,Thermal,bins,bsteps,J2);
figure
plot(S_array(:,2),S_array(:,1),'ro-')

%% 主函数
function [S_array]=mcmc(n,T,Thermal,bins,bsteps,J2)
beta=1/T;
[ising]=generate(n);
[nbor]=neighbor(n);
for j=1:Thermal
    [ising]=flip_spin(ising,nbor,n,beta,J2);
end
[correct]=spatial_association(ising,n);
S_array=zeros(length(correct(:,2)),2);
for j=1:bins
    for i=1:bsteps
        [ising]=flip_spin(ising,nbor,n,beta,J2);
    end
    [correct]=spatial_association(ising,n);
    S_array=S_array+correct;
end
S_array=S_array/bins;
fprintf('temperature is %f\t',T); fprintf('--已完成--\n');
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
%% correction function
function [correct]=spatial_association(spin,n)
array=zeros(n^4-n^2*(1+n^2)/2,2);
count=1;
for i=1:n*n
    iyc= fix((i-1)/n)+1;
    ixc=i-(iyc-1)*n;
    for j=1:(n*n-i)
        iy= fix((j+i-1)/n)+1; ix=j+i-(iy-1)*n;
        array(count,1)=spin(i)*spin(i+j);%关联函数
        x=abs(ixc-ix);y=abs(iyc-iy);
        if x>n/2
            x=n-x;
        end
        if y>n/2
            y=n-y;
        end
        array(count,2)=sqrt(x^2+y^2);%空间距离
        count=count+1;
    end
end
%接下来合并相同距离的项，并求出平均关联
correct(:,2)=unique(array(:,2));
correct(:,1)=0;
for i=1:length(correct)
    list=find(array(:,2)==correct(i,2));
    for j=1:length(list)
        correct(i,1) = correct(i,1)+array(list(j),1);
    end
    correct(i,1)=correct(i,1)/length(list);
end
end