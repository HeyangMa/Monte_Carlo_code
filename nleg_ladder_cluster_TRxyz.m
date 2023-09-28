%metropolis: many-body interaction
%H= \sum Kx*SiSj + Ky*SiSk + Kz*SjSk
% = \sum K_min(SiSj+SiSk+SjSk) + delta_Kx*SiSj + delta_Ky*SiSk + delta_Kz*SjSk    (没用到)
clc
clear
% close all
%% parameter setting
n=6;     %even
T=(0.1:0.1:2);  %修改温度范围
Kx=1.2; Ky=1; Kz=0.8;  %取正为fan铁磁
Thermal=5000;      %弛豫步数
sampling=20000;    %采样数
interval=2;        %采样间隔
%% 预分配
m_sampling=zeros(1,sampling);
e_sampling=zeros(1,sampling);
m=zeros(1,length(T));
m2=zeros(1,length(T));
e=zeros(1,length(T));
cv=zeros(1,length(T));
B_r=zeros(1,length(T));
%% main
[nbor]=neighbor(n);
for t=1:length(T)
    beta=1/T(t);
    B1 = [exp(-beta*Kx), exp(beta*Kx), exp(beta*Kx), exp(-beta*Kx)];
    B1 = reshape(B1, 2, 2);
    B2 = [exp(-beta*Ky), exp(beta*Ky), exp(beta*Ky), exp(-beta*Ky)];
    B2 = reshape(B2, 2, 2);
    B3 = [exp(-beta*Kz), exp(beta*Kz), exp(beta*Kz), exp(-beta*Kz)];
    B3 = reshape(B3, 2, 2);
    [ising]=generate(n);
    
    for j=1:Thermal%弛豫
        [ising]=line(n, nbor, ising, B1,B2,B3);
    end
    for j=1:sampling%采样
        for kk=1:interval
            [ising]=line(n, nbor, ising, B1,B2,B3);
        end
        [e_sampling(j),m_sampling(j)]=calculate(ising,Kx,Ky,Kz,n,nbor);
    end
    %average
    e(t)=mean(e_sampling);
    m(t)=mean(m_sampling);
    cv(t)=beta^2 * var(n^2*e_sampling);
    m2(t)=beta * n^2 * var(m_sampling);
    %宾德累计量
    B_r(t)=mean(m_sampling.^2)/mean(m_sampling);
    
    fprintf('temperature is %f\n',T(t));
end
figure(1);plot(T,e,'ro-');xlabel('T');ylabel('E');hold on
figure(2);plot(T,m,'ro-');xlabel('T');ylabel('M');hold on
figure(3);plot(T,cv,'ro-');xlabel('T');ylabel('C_{v}');hold on
figure(4);plot(T,m2,'ro-');xlabel('T');ylabel('\chi');hold on
figure(5);plot(T,B_r,'ro-');xlabel('T');ylabel('B_{r}');hold on

%% generate ground state ising model
function [ising]=generate(n)
ising=ones(1,n*n);
end
%% generate table of near neighbor
function [nbor]=neighbor(n)
nbor = zeros(n*n,6);
for i=linspace(1,n*n,n*n)
    iy=fix((i-1)/n)+1; %i点的y轴坐标
    ix=i-(iy-1)*n;     %i点的x轴坐标
    ixp = rem(ix+1,n); %右侧点的x坐标
    if (ixp == 0)
        ixp = ixp+n;
    end
    iyp = rem(iy+1,n); %右上侧点的y坐标
    if (iyp == 0)
        iyp = iyp+n;
    end
    ixm = ix-1;         %左侧点的x坐标
    if (ixm == 0)
        ixm=ixm+n;     %周期
    end
    iym = iy-1;         %左下侧点的y坐标
    if (iym == 0)
        iym=iym+n;     %周期
    end
    nbor(i,1) = ixp+(iy-1)*n;  %右
    nbor(i,2) = ix+(iyp-1)*n;  %右上
    nbor(i,3) = ixm+(iyp-1)*n; %左上
    nbor(i,4) = ixm+(iy-1)*n;  %左
    nbor(i,5) = ix+(iym-1)*n;  %左下
    nbor(i,6) = ixp+(iym-1)*n; %右下
end
end
%% flip spin
function [res]=line(n, nbor, spin, B1,B2,B3)

headx = unidrnd(n);
heady = unidrnd(n);
l_site = (heady-1)*n+headx;
tailr = headx+n-2;

%初始化
ispin = 0;
flag1 = -ones(n^2,1);
order = -ones(n^2,1);
stack = -ones(n^2,1);
fm    =  ones(2,n^2);

%堆栈
direct=unidrnd(3);
for i = headx:tailr
    if direct==1
        l_site = nbor(l_site,1);             %点编号
    elseif direct==2
        l_site = nbor(l_site,2);             %点编号
    else
        l_site = nbor(l_site,3);             %点编号
    end
    ispin = ispin + 1;
    stack(ispin) = l_site;                  %存点编号  stack   1            2     3      4      5      6         7              给定一个序号，对应点编号
    flag1(l_site) = 1;                      %                 头的点编号                                      尾的点编号
    order(l_site) = ispin;                  %存线编号  order   头的点编号                                      尾的点编号        给定一个序号，对应线编号
end                                         %                 1            2     3      4      5      6         7
%应注意, 尾巴的线编号是大于头的线编号的,但点编号不一定


for i = ispin:-1:1 %从尾到头遍历
    isk = stack(i);
    temp = [1; 1];
    for j = 1:6  %遍历线上的点的邻居
        
        ik = nbor(isk, j);
        
        if flag1(ik) ~= 1 %非线上的点
            idt = spin(ik);
            if idt == -1
                idt = 2;
            end
            if j==1 || j==4
                I2 = [B1(1, idt); B1(2, idt)];
            elseif j==2 || j==5
                I2 = [B2(1, idt); B2(2, idt)];
            elseif j==3 || j==6
                I2 = [B3(1, idt); B3(2, idt)];
            end
            temp = temp .* I2;
        end
        
        if order(isk) < order(ik) && flag1(ik) > 0  %邻居也是线上的点  并  8-->9
            idt = spin(ik);
            if idt == -1
                idt = 2;
            end
            if j==1 || j==4
                I2 = [B1(1, idt); B1(2, idt)];
            elseif j==2 || j==5
                I2 = [B2(1, idt); B2(2, idt)];
            elseif j==3 || j==6
                I2 = [B3(1, idt); B3(2, idt)];
            end
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
    if direct==1
        I1 = [B1(1, idt); B1(2, idt)];
    elseif direct==2
        I1 = [B2(1, idt); B2(2, idt)];
    elseif direct==3
        I1 = [B3(1, idt); B3(2, idt)];
    end
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
%% calculate_energy
function [e_sampling,m_sampling]=calculate(ising,Kx,Ky,Kz,n,nbor)
energy=0.0;
m_A=0.0;
m_B=0.0;
m_C=0.0;
for j=1:n*n
    energy=energy + Kx * ising(j)*ising(nbor(j,4));
    energy=energy + Ky * ising(j)*ising(nbor(j,2));
    energy=energy + Kz * ising(j)*ising(nbor(j,3));
    iy=fix((j-1)/n)+1; %i点的y轴坐标
    ix=j-(iy-1)*n;     %i点的x轴坐标
    if mod(iy,3)==1 && mod(ix,3)==1 || mod(iy,3)==2 && mod(ix,3)==2 || mod(iy,3)==0 && mod(ix,3)==0
        m_A=m_A+ising(j);
    end
    if mod(iy,3)==1 && mod(ix,3)==0 || mod(iy,3)==2 && mod(ix,3)==1 || mod(iy,3)==0 && mod(ix,3)==2
        m_B=m_B+ising(j);
    end
    if mod(iy,3)==1 && mod(ix,3)==2 || mod(iy,3)==2 && mod(ix,3)==0 || mod(iy,3)==0 && mod(ix,3)==1
        m_C=m_C+ising(j);
    end
end
e_sampling=energy/n^2;
m_sampling=((m_C)^2+(m_B)^2+(m_A)^2)/n^2;
end