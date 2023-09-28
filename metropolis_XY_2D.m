% metropolis algorithm for XY models on the square lattice
% H_{c}=\sum -J*cos(delta_theta)
% author: MHY
% 2022-09
clc
clear
% close all
%% 参数设置
n=5;                    %维度格点数
J=1;                    %耦合常数
T=(0:0.1:2);            %温度设置
Thermal=5000;           %弛豫步数
samples=5000;
bsteps=1;
%% 不同温度
data=zeros(length(T),7);
for i=1:length(T)
    [data(i,:)]=mcmc(n,J,T(i),Thermal,samples,bsteps);
end
%% figure
figure(1);plot(T,data(:,1),'ko-');ylabel('E');  hold on
figure(2);plot(T,data(:,2),'ko-');ylabel('vortex1');  hold on
figure(2);plot(T,data(:,3),'k*-');ylabel('vortex2');  hold on
figure(3);plot(T,data(:,4),'ko-');ylabel('flu vortex1');  hold on
figure(3);plot(T,data(:,5),'k*-');ylabel('flu vortex2');  hold on
figure(4);plot(T,data(:,6),'ko-');ylabel('m1');  hold on
figure(4);plot(T,data(:,7),'k*-');ylabel('m2');  hold on

%% 主函数
function [a1data]=mcmc(n,J,T,Thermal,samples,bsteps)
beta=1.0/T;
[nbor]=neighbor(n);
[xy]=generate_spin(n);
for i=1:Thermal%弛豫
    [xy]=flip_spin(xy,n,nbor,J,beta);
end
a1data=zeros(1,3);
data1=zeros(samples,3);
for i=1:samples%采样
    for j=1:bsteps
        [xy]=flip_spin(xy,n,nbor,J,beta);
    end
    [data1(i,1),data1(i,2),data1(i,3),data1(i,4),data1(i,5)]=calculate(xy, nbor, n, J);
end
a1data(1)=mean(data1(:,1));%energy
a1data(2)=mean(data1(:,2));%vortex1
a1data(3)=mean(data1(:,3));%vortex2
a1data(6)=mean(data1(:,4));%m
a1data(7)=mean(data1(:,5));%m2
a1data(4)=var(data1(:,2));%flu vortex1
a1data(5)=var(data1(:,3));%flu vortex2
%显示当前温度
fprintf('temperature is %f\t',T); fprintf('--已完成--\n');
% draw(xy,n,T,nbor);
end
%% 初始化自旋
function [xy]=generate_spin(n)
xy=rand(n*n,1)*2*pi;
end
%% 邻居关系
function [nbor]=neighbor(n)
nbor=zeros(n*n,4);
for ispin=1:n*n
    iy= fix((ispin-1)/n)+1;
    ix=ispin-(iy-1)*n;
    ixp=ix+1-fix(ix/n)*n;
    iyp=iy+1-fix(iy/n)*n;
    ixm=ix-1+fix((n-ix+1)/n)*n;
    iym=iy-1+fix((n-iy+1)/n)*n;
    nbor(ispin,1)=(iy-1)*n+ixp;%右邻居    1   2  3  4
    nbor(ispin,2)=(iym-1)*n+ix;%上邻居    5   6  7  8
    nbor(ispin,3)=(iy-1)*n+ixm;%左邻居    9  10 11 12
    nbor(ispin,4)=(iyp-1)*n+ix;%下邻居    13 14 15 16
end
end
%% 一个 monte carlo 步
function [xy]=flip_spin(xy,n,nbor,J,beta)
for i=1:n*n
    %rand pick one site from xy
    csite=i;
    right=nbor(csite,1); up=nbor(csite,2);  left=nbor(csite,3);  down=nbor(csite,4);
    %calculate energy (old)
    E_old=-J *cos(xy(csite)-xy(left));              %left
    E_old=E_old - J *cos(xy(csite)-xy(right));      %right
    E_old=E_old - J *cos(xy(csite)-xy(up));         %up
    E_old=E_old - J *cos(xy(csite)-xy(down));       %down
    %attempt flip spin
    alpha=rand*2*pi;
    angle =2*alpha-xy(csite);
    while angle>2*pi || angle<0
        if angle>2*pi
            angle=angle-2*pi;
        end
        if angle<0
            angle=angle+2*pi;
        end
    end
    %calculate energy (new)
    E_new=-J *cos(angle-xy(left));               %left
    E_new=E_new - J *cos(angle-xy(right));       %right
    E_new=E_new - J *cos(angle-xy(up));          %up
    E_new=E_new - J *cos(angle-xy(down));        %down
    %calculate delta energy
    deltaE=E_new-E_old;
    if rand<exp(-1*deltaE*beta)
        xy(csite)=angle;
    end
end
end
%% 计算物理量
function [energy,vortex1,vortex2,m,m2]=calculate(xy, nbor, n, J)
energy = 0.0;
for j = 1:n*n
    right=nbor(j,1); up=nbor(j,2);
    energy=energy - J *cos(xy(j)-xy(right));           %right
    energy=energy - J *cos(xy(j)-xy(up));              %up
end
energy = energy / (n * n);
vortex1=count_integer_vortex(xy, nbor,n);
xy2=xy*2;
for i = 1:n^2
    while xy2(i) > 2 * pi || xy2(i) < 0
        if xy2(i) > 2 * pi
            xy2(i) = xy2(i) - 2 * pi;
        end
        if xy2(i) < 0
            xy2(i) = xy2(i) + 2 * pi;
        end
    end
end
vortex2=count_integer_vortex(xy2, nbor,n);
m=sqrt(sum(cos(xy))^2 + sum(sin(xy))^2)/n^2;
m2=sqrt(sum(cos(xy/2))^2 + sum(sin(xy/2))^2)/n^2;
end
%% calculate vortex
function [vortex]=count_integer_vortex(xy, nbor,n)
%---------integer vortex---------
theta = zeros( 4, 1);
v = zeros( 4, 1);
vortex = 0;
for i = 1:n^2
    theta(1) = xy(nbor(nbor(i, 1), 2));
    theta(2) = xy(nbor(i, 2));
    theta(3) = xy(i);
    theta(4) = xy(nbor(i, 1));
    v(1) = theta(2) - theta(1);
    v(2) = theta(3) - theta(2);
    v(3) = theta(4) - theta(3);
    v(4) = theta(1) - theta(4);
    for k = 1:4
        if v(k) <= -pi
            v(k) = v(k) + 2 * pi;
        end
        if v(k) > pi
            v(k) = v(k) - 2 * pi;
        end
    end
    kv = sum(v) / 2 / pi;
    if  floor(kv)==1
        vortex = vortex + 1;
    end
    if floor(kv)==-1
        vortex = vortex + 1;
    end
end
end
%% 画涡旋
function draw(xy,n,T,nbor)
%----------  1 倍 theta  ----------
sin_xy=sin(xy);
cos_xy=cos(xy);
figure;hold on
for i=1:n*n
    iy=fix((i-1)/n)+1; ix=i-(iy-1)*n;
    quiver(ix-cos_xy(i)/3,iy-sin_xy(i)/3,cos_xy(i)/3*2,sin_xy(i)/3*2,'MaxHeadSize',2,'AutoScaleFactor',0.5,'AutoScale','off','color',[174,174,174]/255);
    %     text(ix,iy+0.1,num2str(i));text(ix,iy+0.3,num2str(xy(i)));
end
grid on;
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);  set(gca,'ticklength',[0 0]);
canvasX=2;      %画布右下角点的X坐标
canvasY=2;      %画布右下角点的Y坐标
canvasL=18;     %画布宽
canvasH=15;    %画布高
set(gcf,'unit','centimeters','position',[canvasX canvasY canvasL canvasH]) %定义值的单位为厘米
title(['T=',num2str(T)]);
error=0.2;
for i=1:n*n
    iy=fix((i-1)/n)+1; ix=i-(iy-1)*n;
    theta(1)=xy(i);
    theta(2)=xy(nbor(i,1));        %right           4       3
    theta(3)=xy(nbor(nbor(i,1),4));%right up
    theta(4)=xy(nbor(i,4));        %up              1       2
    v(1)=theta(2)-theta(1);
    v(2)=theta(3)-theta(2);
    v(3)=theta(4)-theta(3);
    v(4)=theta(1)-theta(4);
    for k=1:4
        if v(k)>=pi
            v(k)=v(k)-2*pi;
        end
        if v(k)<=-pi
            v(k)=v(k)+2*pi;
        end
    end
    kv=sum(v)/2/pi;
    if  (1-error) < kv && kv < (1+error)  % 1-0.001 < kv < 1+0.001
        text(ix+0.5,iy+0.5,'+','color','b','fontsize',20)
    end
    if (-1-error) < kv && kv < (-1+error) % -1-0.001 < kv < -1+0.001
        text(ix+0.5,iy+0.5,'-','color','m','fontsize',20)
    end
%     if  floor(kv)==1
%         text(ix+0.5,iy+0.5,'+','color','b','fontsize',20)
%     end
%     if floor(kv)==-1
%         text(ix+0.5,iy+0.5,'-','color','m','fontsize',20)
%     end
end
end