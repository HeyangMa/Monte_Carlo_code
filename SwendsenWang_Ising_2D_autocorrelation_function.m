%SW method for ising model: many-body interaction
%H=\sum -J2*SiSj - J4*SiSjSmSn - J6*SiSjSmSnSpSq
%T=(2/log(1+sqrt(2)));
%% 清理内存
clc
clear
% close all
%% 参数设置
n=8;
T=(0:0.2:5);
J2=1;J4=0.2;J6=0;
Thermal=1000;          %弛豫步数
bins=100;             %bins数目
bsteps=1000;           %每个bin的步数
ntime=50;              %算相变点处的自关联函数, 显示的时间( 1 2 3 4 5 ... 100)
%% 不同温度模拟
for i=1:length(T)
    [e,m,cv,ms,e_Variance,m_Variance,~,~]=mcmc(n,T(i),Thermal,bins,bsteps,J2,J4,J6);
    figure(1);
    subplot(2,2,1); errorbar(T(i),e,e_Variance,'ro:');  xlabel('T'); ylabel('E');hold on
    subplot(2,2,2); errorbar(T(i),m,m_Variance,'ro:');  xlabel('T'); ylabel('M');hold on
    subplot(2,2,3); plot(T(i),cv,'ro');  xlabel('T'); ylabel('C_{v}');hold on
    subplot(2,2,4); plot(T(i),ms,'ro');  xlabel('T'); ylabel('\chi'); hold on
end
%% 计算自关联函数
% T=(2/log(1+sqrt(2)));
% T=2.4955;
% [~,~,~,~,~,~,e_steps,m_steps]=mcmc(n,T,Thermal,bins,bsteps,J2,J4,J6);
% [AC_function,AC_errorbar]=auto_correlation_function(m_steps,ntime,bins,bsteps);
% figure(2)
% errorbar(0:(ntime-1),AC_function,AC_errorbar,'ko-');  xlabel('t'); ylabel('A_{m}(t)');
% xlim([0 50])
% set(gca, 'YScale', 'log')
% ylim([1.5*10^(-3) 1])
% hold on

% data=[T',e',m',cv',ms',Variance'];
% save C:\Users\LENOVO\Desktop\e6.dat -ascii  data
% h=figure(1);
% saveas(h,'C:\Users\LENOVO\Desktop\figure_6.pdf');


%% 主函数
function [e,m,cv,ms,e_Variance,m_Variance,e_steps,m_steps]=mcmc(n,T,Thermal,bins,bsteps,J2,J4,J6)
beta=1/T;
nbor=neighbor(n);%近邻
ising=generate_ising(n);%生成
fprintf('temperature is %f\t',T); fprintf('--已开始--\n');
for j=1:Thermal%弛豫
    [bonds]=configure_bonds(ising,nbor,n,J2,J4,J6,beta);
    [cluster]=search_cluster(bonds,nbor,n);%print_cluster(cluster,n,ising);检查集团并画图
    [ising]=flip_spin(cluster,ising);%一步mc
end
for j=1:bins*bsteps%采样
    [bonds]=configure_bonds(ising,nbor,n,J2,J4,J6,beta);
    [cluster]=search_cluster(bonds,nbor,n);%[cluster_max(j)]=find_max_cluster(cluster);%找出最多的元素数量
    [ising]=flip_spin(cluster,ising);%一步mc
    [e_steps(j),m_steps(j)]=calculate_energy(ising,n,nbor,J2,J4,J6);
end
%mean_cluster_max(t)=mean(cluster_max);
e=mean(e_steps);e_Variance=error_bar(e_steps,bins);
m=mean(m_steps);m_Variance=error_bar(m_steps,bins);
[cv,ms]=calculate_2(e_steps,m_steps,bins,bsteps,beta);
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
%% 标准差
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
%% 画出cluster图
function print_cluster(cluster,n,ising)
cluster=reshape(cluster,n,n)';
%输出为图像
imagesc(cluster)
hold on
for ii=1:n*n
    iy=fix((ii-1)/n)+1;
    ix=ii-(iy-1)*n;
    text(ix,iy,num2str(ising(ii)));
    hold on
end
grid on;set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
for  x=0.5:1:n
    plot([x,x],[0,n+1],'k-')
end
for  y=0.5:1:n
    plot([0,n+1],[y,y],'k-')
end

set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'ticklength',[0 0])
axis square
end
%% 找出最多的元素数量
function [cluster_max]=find_max_cluster(cluster)
cluster=cluster(:);
table = tabulate(cluster);
[~,idx] = max(table(:,2));%[maxCount,idx] = max(table(:,2));
%获取出现次数最多的元素
cluster_max=sum(cluster(:)==table(idx));
end
%% generate fe ising model
function [ising]=generate_ising(n)
ising=ones(1,n*n);
end
%% generate table of near neighbor
function [nbor]=neighbor(n)
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
%% configure bonds
function [bonds]=configure_bonds(ising,nbor,n,J2,J4,J6,beta)
bonds=zeros(n*n,10);
%----------------
for i=1:n*n
    pro_bonds_right = 1 - exp( -2* J2 * beta * ising(i) * ising(nbor(i,1)) );
    pro_bonds_up    = 1 - exp( -2* J2 * beta * ising(i) * ising(nbor(i,2)) );
    if rand<pro_bonds_right
        bonds(i,1)=1;
    end
    if rand<pro_bonds_up
        bonds(i,2)=1;
    end
end
%----------------
for j=1:n*n
    product=ising(j) * ising(nbor(j,1)) * ising(nbor(j,2)) * ising(nbor(nbor(j,2),1));
    pro_bonds_square = 1 - exp( -2* J4 * beta * product );
    
    bondsum=bonds(j,1) + bonds(j,2) + bonds(nbor(j,2),1) + bonds(nbor(j,1),2);
    
    if J4>0
        select=1;
    else
        select=-1;
    end
    
    if product==select
        if bondsum==0
            unid=unidrnd(3);
            if unid==1  && rand<pro_bonds_square
                bonds(j,1)=1;                    %    -----
                bonds(nbor(j,2),1)=1;            %    -----
            elseif unid==2 && rand<pro_bonds_square
                bonds(j,2)=1;                    %    |   |
                bonds(nbor(j,1),2)=1;            %    |   |
            elseif unid==3 && rand<pro_bonds_square
                bonds(j,3)=1;                    %    \  /
                bonds(nbor(j,1),4)=1;            %    /  \
            end
        elseif bondsum==1
            if bonds(j,1)==1  && rand<pro_bonds_square
                bonds(nbor(j,2),1)=1;
            elseif bonds(nbor(j,2),1)==1  && rand<pro_bonds_square
                bonds(j,1)=1;
            elseif bonds(j,2)==1  && rand<pro_bonds_square
                bonds(nbor(j,1),2)=1;
            elseif bonds(nbor(j,1),2)==1  && rand<pro_bonds_square
                bonds(j,2)=1;
            end
        elseif bondsum==2 && rand<pro_bonds_square
            bonds(j,1)=1;
            bonds(j,2)=1;
            bonds(nbor(j,2),1)=1;
            bonds(nbor(j,1),2)=1;
        end
    end
end
%----------------
for i=1:n*n
    
    site_right=nbor(i,1);                       %右
    site_up=nbor(i,2);                          %上
    site_rightup=nbor(site_up,1);               %右上
    site_upup=nbor(site_up,2);                  %上上
    site_upupright=nbor(site_upup,1);           %上上右
    
    product=ising(i)*ising(site_right)*ising(site_up)*ising(site_rightup)*ising(site_upup)*ising(site_upupright);
    
    pro_bonds_J6= 1 - exp( -2* J6 * beta * product );
    
    %|------|
    %  \  /
    %  /  \
    %|------|
    %  \  /
    %  /  \
    %|------|
    
    %3x2
    if product==1  && rand<pro_bonds_J6
        switch 3 %三种放棒方式
            case 1
                
                %在n=6的倍数时使用
                iy= fix((i-1)/n)+1;
                if mod(iy,3)==0
                    bonds(i,2)=1;  bonds(nbor(i,1),2)=1;  bonds(nbor(nbor(i,2),2),1)=1;
                elseif mod(iy,3)==2
                    bonds(i,3)=1;  bonds(nbor(i,2),1)=1;  bonds(nbor(i,1),3)=1;
                elseif mod(iy,3)==1
                    bonds(i,1)=1;  bonds(nbor(i,2),2)=1;  bonds(nbor(nbor(i,1),2),2)=1;
                end
                
            case 2
                
                bonds(cen,1)=1;
                bonds(site_up,1)=1;
                bonds(site_upup,1)=1;
                
            case 3
                
                switch unidrnd(2)
                    case 1
                        bonds(i,1)=1; bonds(nbor(i,2),2)=1; bonds(nbor(nbor(i,2),1),2)=1;
                    case 2
                        bonds(i,2)=1; bonds(nbor(i,1),2)=1; bonds(nbor(nbor(i,2),2),1)=1;
                end
                
        end
        
    end
    
end
end
%% apply the H-K algorithm and search cluster with PBC
function [cluster]=search_cluster(bonds,nbor,n)
cluster=zeros(1,n*n);
largest_label=0;
label=(1:n*n);
for j=1:n*n
    %1:左  2:上  3:左上  4:右上
    %5:左左      6:上上
    %7:左上上    8:右上上      9:左左上     10:右右上
    cen=j;
    xb=bonds(j,1);
    yb=bonds(j,2);
    lb=bonds(j,3);
    rb=bonds(j,4);
    llb=bonds(j,5);
    uub=bonds(j,6);
    luub=bonds(j,7);
    ruub=bonds(j,8);
    llub=bonds(j,9);
    rrub=bonds(j,10);
    bondlist=[xb,yb,lb,rb,llb,uub,luub,ruub,llub,rrub];
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
[cluster]=relabel(n*n,cluster);
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
%search nbor site
function sitenbor=searchnbor(number,nbor,cen)
if number==1
    sitenbor=nbor(cen,1);
elseif number==2
    sitenbor=nbor(cen,2);
elseif number==3
    sitenbor=nbor(nbor(cen,2),1);
elseif number==4
    sitenbor=nbor(nbor(cen,2),3);
    % elseif number==5
    %     sitenbor=nbor(nbor(nbor(cen,3),3),2);
    % elseif number==6
    %     sitenbor=nbor(nbor(nbor(cen,3),2),2);
    % elseif number==7
    %     sitenbor=nbor(nbor(nbor(cen,1),2),2);
    % elseif number==8
    %     sitenbor=nbor(nbor(nbor(cen,1),1),2);
    % elseif number==9
    %     sitenbor=nbor(nbor(cen,2),2);
    % elseif number==10
    %     sitenbor=nbor(nbor(cen,3),3);
end
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
%% flip cluster's spins
function [ising]=flip_spin(cluster,ising)
for i=1:max(max(cluster))
    if rand>0.5 %成功翻转
        [n]=find(cluster==i);
        for j=1:length(n)
            ising(n(j))=-ising(n(j));
        end
    end
end
end
%% calculate
function [energy,m]=calculate_energy(ising,n,nbor,J2,J4,J6)
energy=0.0;
for i=1:n*n
    right=nbor(i,1);up=nbor(i,2);rightup=nbor(right,2);
    upup=nbor(up,2);rightupup=nbor(rightup,2);
    energy=energy - J2 * ising(i) * ( ising(up) + ising(right) );
    energy=energy - J4 * ising(i) * ising(up) * ising(rightup) * ising(right);
    energy=energy - J6 * ( ising(i) * ising(up)   * ising(upup) * ising(rightupup) * ising(rightup) * ising(right) );
end
energy =energy/n^2;
m=abs(sum(sum(ising)))/(n*n);
end
%% 计算比热、磁化率
function [cv,ms]=calculate_2(e,m,bins,bsteps,beta)
j=1;
e2=e.*e;
m2=m.*m;
e_bins=zeros(1,bins);
m_bins=zeros(1,bins);
e2_bins=zeros(1,bins);
m2_bins=zeros(1,bins);
for i=1:length(e)
    e_bins(j)=e_bins(j)+e(i);e2_bins(j)=e2_bins(j)+e2(i);
    m_bins(j)=m_bins(j)+m(i);m2_bins(j)=m2_bins(j)+m2(i);
    while i==bsteps*j
        e_bins(j)=e_bins(j)/bsteps;e2_bins(j)=e2_bins(j)/bsteps;
        m_bins(j)=m_bins(j)/bsteps;m2_bins(j)=m2_bins(j)/bsteps;
        j=j+1;
    end
end
% e_bins=e_bins*n*n;
% m_bins=m_bins*n*n;

cv=beta^2 * (mean(e2_bins)-mean(e_bins)^2);
ms=beta   * (mean(m2_bins)-mean(m_bins)^2);
end