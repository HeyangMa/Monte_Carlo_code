%metropolis: many-body interaction
%H=\sum -J2*SiSj
%T=(2/log(1+sqrt(2)));
%% 清理内存
clc
clear
% close all
%% 参数设置
n=10;
T=(0.1:0.1:4);
J2=1;
Thermal=1000;          %弛豫步数
bins=100;              %bins数目
bsteps=100;            %每个bin的步数
%% 主函数
[nbor]=neighbor(n);
for t=1:length(T)
    beta=1/T(t);
    B = [exp(beta), exp(-beta), exp(-beta), exp(beta)];     %  e^(-beta<-J>)      e^(-beta<J>)
    B = reshape(B, 2, 2); 
    [ising]=generate(n);                                    %  e^(-beta<J>)      e^(-beta<-J>)
    for j=1:Thermal
        [ising]=line(n, nbor, ising, B);
    end
    for j=1:bins*bsteps
            [ising]=line(n, nbor, ising, B);
            [e(j)]=calculate(ising,nbor,n,J2);
            m(j)=abs(sum(ising))/n;
    end
    energy(t)=mean(e);
    mag(t)=mean(m);
    fprintf('temperature is %f\t',T(t)); fprintf('--已完成--\n');
end
figure(1);hold on;plot(T,energy,'g>-');
figure(2);hold on;plot(T,mag,'g>-');
%% 初始构型
function [ising]=generate(n)
ising=ones(n,1);
end
%% 邻居表
function [nbor]=neighbor(n)
nbor=zeros(n,2);
for ispin=1:n
    nbor(ispin,1)=ispin+1;%右邻居
    nbor(ispin,2)=ispin-1;%左邻居
    if ispin==1
        nbor(ispin,2)=n;  %左邻居
    end
    if ispin==n
        nbor(ispin,1)=1;  %右邻居
    end
end
end
%% 一步mc
function [res]=line(n, nbor, spin, B)

headx = unidrnd(n);

tailr = headx+n-2;

%初始化
ispin = 0;
flag1 = -ones(n,1);
order = -ones(n,1);
stack = -ones(n,1);
fm    =  ones(2,n);

%堆栈
for i = headx:tailr
    linex = mod(i, n);
    if linex==0
        linex=n;
    end
    l_site = linex;                         %点编号
    ispin = ispin + 1;
    stack(ispin) = l_site;                  %存点编号  stack   1            2     3      4      5      6         7              给定一个序号，对应点编号
    flag1(l_site) = 1;                      %                 头的点编号                                      尾的点编号
    order(l_site) = ispin;                  %存线编号  order   头的点编号                                      尾的点编号        给定一个序号，对应线编号
end                                         %                 1            2     3      4      5      6         7  
%应注意, 尾巴的线编号是大于头的线编号的,但点编号不一定


for i = ispin:-1:1 %从尾到头遍历
    isk = stack(i);
    temp = [1; 1];
    for j = 1:2  %遍历线上的点的邻居
        
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
function [energy]=calculate(ising,nbor,n,J2)
energy=0.0;

for i=1:n
    
    energy=energy - J2 * ising(i) * ( ising(nbor(i,1))  );

end

energy=energy/(n);
end