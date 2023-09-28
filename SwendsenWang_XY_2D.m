%SW method for ising-xy model: four-body interaction
%H=\sum -J2 * Si * Sj * \vec{Si} * \vec{Sj}
% close all
clear
clc
%% parameter setting
n=6;
T=(0:0.1:2);
J2=1;
Thermal=5000;   %弛豫步数
bins=5000;
bsteps=1;
%% 预分配
e_samp = zeros(1,bins*bsteps);
m_samp = zeros(1,bins*bsteps);
e   = zeros(1,length(T));
cv  = zeros(1,length(T));
m   = zeros(1,length(T));
chi = zeros(1,length(T));
R2  = zeros(1,length(T));
e_error  = zeros(1,length(T));
m_error  = zeros(1,length(T));
cv_error = zeros(1,length(T));
ms_error = zeros(1,length(T));
R2_error = zeros(1,length(T));
%% main
nbor=neighbor(n);%近邻
spin_initial=generate_spin(n);

for t=1:length(T)
    beta=1/T(t);
    spin=spin_initial;
    for i=1:Thermal
        angle=rand*2*pi;
        [bonds]=configure_bonds(spin,nbor,n,J2,angle,beta);
        [cluster]=search_cluster(bonds,nbor,n);
        [spin]=flip_spin(cluster,spin,angle);
    end
    for i=1:bsteps*bins
        angle=rand*2*pi;
        [bonds]=configure_bonds(spin,nbor,n,J2,angle,beta);
        [cluster]=search_cluster(bonds,nbor,n);
        [spin]=flip_spin(cluster,spin,angle);
        [e_samp(i),m_samp(i)]=calculate(spin,nbor,n,J2);
    end
    %energy
    e(t)=mean(e_samp);
    e_error(t)=error_bar(e_samp,bins);
    %mag
    m(t)=mean(m_samp);
    m_error(t)=error_bar(m_samp,bins);
    %cv
    cv(t) =beta^2 * var(n^2*e_samp);
    %ms
    chi(t)=beta * n^2 * var(m_samp);
    %binder
    R2(t)=mean(m_samp.^2)^2/mean(m_samp.^4);
    %error
    [cv_error(t), ms_error(t), R2_error(t)]=other_errorbar(e_samp, m_samp, n, T(t), bins);
    fprintf('temperature is %f\n',T(t));
end
%%
% figure(1);hold on;errorbar(T,e,e_error,'ko:');    xlabel('T'); ylabel('E');
% figure(2);hold on;errorbar(T,m,m_error,'ko:');    xlabel('T'); ylabel('M');
% figure(3);hold on;errorbar(T,cv,cv_error,'ko:');  xlabel('T'); ylabel('C_{v}');
% figure(4);hold on;errorbar(T,chi,ms_error,'ko:'); xlabel('T'); ylabel('\chi');
% figure(5);hold on;errorbar(T,R2,R2_error,'ko:');  xlabel('T'); ylabel('R2');
%%
figure(1);hold on;plot(T,e,  'ko:');  xlabel('T'); ylabel('E');
% figure(2);hold on;plot(T,m,  'ko:');  xlabel('T'); ylabel('M');
% figure(3);hold on;plot(T,cv, 'ko:');  xlabel('T'); ylabel('C_{v}');
% figure(4);hold on;plot(T,chi,'ko:');  xlabel('T'); ylabel('\chi');
% figure(5);hold on;plot(T,R2, 'ko:');  xlabel('T'); ylabel('R2');

%% generate ising
function [spin]=generate_spin(n)
spin=zeros(n*n,1);
end
%% neighbor
function [nbor]=neighbor(n)
nbor=zeros(n*n,4);
for ispin=1:n*n
    iy=fix((ispin-1)/n)+1;
    ix=ispin-(iy-1)*n;
    ixp=ix+1-fix(ix/n)*n;
    iyp=iy+1-fix(iy/n)*n;
    ixm=ix-1+fix((n-ix+1)/n)*n;
    iym=iy-1+fix((n-iy+1)/n)*n;
    nbor(ispin,1)=(iy-1)*n+ixp;%右邻居
    nbor(ispin,2)=(iyp-1)*n+ix;%上邻居
    nbor(ispin,3)=(iy-1)*n+ixm;%左邻居
    nbor(ispin,4)=(iym-1)*n+ix;%下邻居
end
end
%% [c,d]=project_t(a,b)
function [proj_value,d]=project_t(a,b)
% 给定角度和投影的目标角度
theta = a;
theta_proj = b;
% 计算角向量和投影方向向量
vec_angle = [cos(theta), sin(theta)];
proj_dir = [cos(theta_proj), sin(theta_proj)];
% 计算投影值
proj_value = dot(vec_angle, proj_dir) / norm(proj_dir);
% 给定两个角度和方向
theta1 = a;
theta2 = b;
dir_len = 1; %方向向量长度

% 计算两个方向向量
dir1 = [dir_len * cos(theta1), dir_len * sin(theta1)];
dir2 = [dir_len * cos(theta2), dir_len * sin(theta2)];

% 计算反射向量
reflection_dir = dir1 - (2 * dot(dir1, dir2) / (norm(dir2)^2)) * dir2;

% 计算反射角
theta_refl = atan2(reflection_dir(2), reflection_dir(1));
if theta_refl < 0
    theta_refl = 2 * pi + theta_refl; % 保证反射角度在0到360之间
end
angle_refl = theta_refl * 180 / pi;
d=angle_refl * pi / 180;
end
%% configure bonds
function [bonds]=configure_bonds(spin,nbor,n,J2,angle,beta)
bonds=zeros(n*n,2);
%计算投影
alpha=zeros(n^2,1);
for j=1:n*n
    [alpha(j),~]=project_t(spin(j),angle);
end
%投影后按ising自旋配置bonds
for j=1:n*n
    %1
    pro_bonds = max([0, 1 - exp( -2 * J2 * beta * alpha(j) * alpha(nbor(j,1)) )]);
    if rand<pro_bonds
        bonds(j,1)=1;
    end
    %2
    pro_bonds = max([0, 1 - exp( -2 * J2 * beta * alpha(j) * alpha(nbor(j,2)) )]);
    if rand<pro_bonds
        bonds(j,2)=1;
    end
end
end
%% apply the H-K algorithm and search cluster with PBC
function [cluster]=search_cluster(bonds,nbor,n)
cluster=zeros(1,n*n);
largest_label=0;
label=(1:n*n);
for j=1:n*n
    cen=j;
    bondlist=[bonds(j,1),bonds(j,2)];
    bondsum=sum(bondlist);
    if cluster(j)==0%-------------------------------------当前点无标签
        if bondsum==0               %无bond
            largest_label=largest_label+1;
            cluster(j)=largest_label;
        elseif bondsum==1           %一个bond
            number=find(bondlist==1);
            sitenbor=searchnbor(number,nbor,cen);
            if cluster(sitenbor)==0
                largest_label=largest_label+1;
                cluster(j)=largest_label;
                cluster(sitenbor)=largest_label;
            else
                [cluster(j),label]=find_root(cluster(sitenbor),label);
            end
        elseif bondsum>1            %多个bond
            numberlist=find(bondlist==1);
            varlist=zeros(1,length(numberlist));
            sitenbor1=zeros(1,length(numberlist));
            for k=1:length(numberlist)
                sitenbor1(k)=searchnbor(numberlist(k),nbor,cen);
                varlist(k)=cluster(sitenbor1(k)); %存贮cluster标签
            end
            labellist=find(varlist~=0);
            if sum(labellist)==0
                largest_label=largest_label+1;
                for k=1:length(numberlist)
                    cluster(sitenbor1(k))=largest_label;
                end
                cluster(j)=largest_label;
            else
                notzerolabellist=varlist(varlist~=0);
                [label,minlabel]=union(notzerolabellist,label);
                for k=1:length(numberlist)
                    cluster(sitenbor1(k))=minlabel;
                end
                cluster(j)=minlabel;
            end
        end
    else%-----------------------------------------------当前点有标签
        if bondsum==0               %无bond
            continue
        elseif bondsum>0          %多个bond(自身带一个)
            numberlist=find(bondlist==1);
            varlist=zeros(1,length(numberlist));
            sitenbor1=zeros(1,length(numberlist));
            for k=1:length(numberlist)
                sitenbor1(k)=searchnbor(numberlist(k),nbor,cen);
                varlist(k)=cluster(sitenbor1(k)); %存贮cluster标签
            end
            labellist=find(varlist~=0);
            if sum(labellist)==0
                for k=1:length(numberlist)
                    [a,label]=find_root(cluster(j),label);
                    cluster(sitenbor1(k))=a;
                end
            else
                notzerolabellist=varlist(varlist~=0);
                [label,minlabel]=union(notzerolabellist,label);
                [a,label]=find_root(cluster(j),label);
                sminlabel=min(minlabel,a);
                label(minlabel)=sminlabel;    label(a)=sminlabel;
                for k=1:length(numberlist)
                    cluster(sitenbor1(k))=minlabel;
                end
                cluster(j)=minlabel;
            end
        end
    end
end
for j=1:n*n
    [cluster(j),label]=find_root(cluster(j),label);
end
[cluster]=relabel(n,cluster);
end
%search nbor site
function sitenbor=searchnbor(number,nbor,cen)
if number==1
    sitenbor=nbor(cen,1);
elseif number==2
    sitenbor=nbor(cen,2);
    % elseif number==3
    %     sitenbor=nbor(cen,5);
    % elseif number==4
    %     sitenbor=nbor(nbor(cen,5),1);
    % elseif number==5
    %     sitenbor=nbor(nbor(cen,5),3);
    % elseif number==6
    %     sitenbor=nbor(nbor(cen,5),2);
    % elseif number==7
    %     sitenbor=nbor(nbor(cen,5),4);
    % elseif number==8
    %     sitenbor=nbor(nbor(nbor(cen,1),1),2);
    % elseif number==9
    %     sitenbor=nbor(nbor(cen,2),2);
    % elseif number==10
    %     sitenbor=nbor(nbor(cen,3),3);
end
end
%find root label
function [y,label]=find_root(m,label)
y=m;
while label(y)~=y
    y=label(y);
end
end
%union bonds
function [label,m]=union(notzerolabellist,label)
var=zeros(1,length(notzerolabellist));
for k=1:length(notzerolabellist)
    [var(k),label]=find_root(notzerolabellist(k),label);
end
a=min(var);
for k=1:length(notzerolabellist)
    label(var(k))=a;
end
m=a;
end
%relabel cluster
function [cluster]=relabel(n,cluster)
relabel_number=1;
for i=1:n*n
    if cluster(cluster==i)~=0
        cluster(cluster==i)=relabel_number;
        relabel_number=relabel_number+1;
    end
end
end
%% cluster flip spins
function [spin]=flip_spin(cluster,spin,angle)
%flip
for i=1:max(max(cluster))
    if rand>0.5 %成功翻转
        [n]=find(cluster==i);
        for j=1:length(n)
            %xy
            [~,spin(n(j))]=project_t(spin(n(j)),angle);
        end
    end
end
end
%% calculate
function [energy,m]=calculate(spin,nbor,n,J2)
energy=0.0;
for i=1:n*n
    energy = energy - J2 * cos(spin(i) - spin(nbor(i,1)));
    energy = energy - J2 * cos(spin(i) - spin(nbor(i,2)));
end
energy=energy/(n*n);
m=sum(cos(spin))/n^2;
end
%% calculate other standard deviation
function [cv_error, ms_error, R2_error]=other_errorbar(e, m, n, T, bins)
cv_bin = zeros(bins, 1);
ms_bin = zeros(bins, 1);
R2_bin = zeros(bins, 1);
for i = 1:bins
    e_bin = e(1+(i-1)*(length(e) / bins):i*(length(e) / bins));
    m_bin = m(1+(i-1)*(length(m) / bins):i*(length(m) / bins));
    %calculate specific heat
    cv_bin(i) = var(e_bin) / T^2;
    %calculate magnetic susceptibility
    ms_bin(i) = n^2 * (sum(m_bin .^ 2) / length(m_bin) - (sum(m_bin) / length(m_bin))^2) / T;
    %calculate binder  cumulant
    m_4 = sum(m_bin .^ 4) / length(m_bin);
    m_2 = sum(m_bin .^ 2) / length(m_bin);
    R2_bin(i) = m_4 / m_2^2;
end
cv_error = std(cv_bin);
ms_error = std(ms_bin);
R2_error = std(R2_bin);
end
%% error_bar
function y=error_bar(x,bins)
j=1;
bin=zeros(1,length(x)/bins);
for i=1:length(x)
    bin(j)=bin(j)+x(i);
    while i==bins*j
        bin(j)=bin(j)/bins;
        j=j+1;
    end
end
y=0;
for i=1:length(bin)
    y=y+(bin(i)-mean(bin))^2;
end
y=sqrt(y/(length(bin)*(length(bin)-1)));
end