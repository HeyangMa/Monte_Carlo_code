%metropolis: many-body interaction
%H=\sum -J2*SiSj
%T=(2/log(1+sqrt(2)));
%% 清理内存
clc;clear;
% close all
%% 参数设置
n=4;
T=(0.1:0.1:8);
J2=1;
Thermal=1000;          %弛豫步数
bins=100;              %bins数目
bsteps=1000;           %每个bin的步数
%% 主函数
for t=1:length(T)
    beta=1/T(t);
    B = [exp(beta), exp(-beta), exp(-beta), exp(beta)];
    B = reshape(B, 2, 2);
    [ising]=generate(n);
    [nbor]=neighbor(n);
    for j=1:Thermal
        [ising]=line(n, nbor, ising, B);
    end
    for j=1:bins*bsteps
            [ising]=line(n, nbor, ising, B);
            [e(j),m(j)]=calculate(ising,nbor,n,J2);
    end
    energy(t)=mean(e);
    mag(t)=mean(m);
    fprintf('temperature is %f\t',T(t)); fprintf('--已完成--\n');
end
figure(1);hold on;plot(T,energy,'b*-');
figure(2);hold on;plot(T,mag,'b*-');
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
function [res]=line(n, nbor, spin, B)

headx = unidrnd(n);
heady = unidrnd(n);

tailr = headx+n-2;

%初始化
ispin = 0;
flag1 = -ones(n^2,1);
order = -ones(n^2,1);
stack = -ones(n^2,1);
fm    =  ones(2,n^2);

%堆栈
for i = headx:tailr
    linex = mod(i, n);
    if linex==0
        linex=n;
    end
    l_site = (heady-1)*n+linex;             %点编号
    ispin = ispin + 1;
    stack(ispin) = l_site;                  %存点编号  stack   1            2     3      4      5      6         7              给定一个序号，对应点编号
    flag1(l_site) = 1;                      %                 头的点编号                                      尾的点编号
    order(l_site) = ispin;                  %存线编号  order   头的点编号                                      尾的点编号        给定一个序号，对应线编号
end                                         %                 1            2     3      4      5      6         7  
%应注意, 尾巴的线编号是大于头的线编号的,但点编号不一定


for i = ispin:-1:1 %从尾到头遍历
    isk = stack(i);
    temp = [1; 1];
    for j = 1:4  %遍历线上的点的邻居
        
        ik = nbor(isk, j);
        
        if flag1(ik) ~= 1 %非线上的点
            idt = spin(ik);
            if idt == -1
                idt = 2;
            end
            I2 = [B(1, idt); B(2, idt)];
            temp = temp .* I2;
        end
        
        if order(isk) < order(ik) && flag1(ik) > 0  %邻居也是线上的点  并  8-->9
            idt = spin(ik);
            if idt == -1
                idt = 2;
            end
            I2 = [B(1, idt); B(2, idt)];
            temp = temp .* I2;
        end
    end
    fm(:, stack(i)) = temp;
end

%采样
i = 1;
ispin = stack(i);
pb = fm(:,ispin);
if rand < (pb(1) / sum(pb))
    spin(ispin) = 1;
else
    spin(ispin) = -1;
end

for i = 2:n-1
    idt = spin(stack(i-1));
    if idt == -1
        idt = 2;
    end
    I1 = [B(1, idt); B(2, idt)];
    temp = fm(:, stack(i));
    pb = I1 .* temp;
    
    ispin = stack(i);
    
    if rand < (pb(1) / sum(pb))
        spin(ispin) = 1;
    else
        spin(ispin) = -1;
    end
end
res = spin;
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