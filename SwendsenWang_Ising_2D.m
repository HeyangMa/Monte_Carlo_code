%SW method for ising model: many-body interaction
%H=\sum -J2*SiSj - J4*SiSjSmSn - J6*SiSjSmSnSpSq
%--------------------------------------------------------------------------
%parameter setting
function SwendsenWang_Ising_2D
n=6;
Nsite=n*n;
T=(0:0.2:4);
beta=1./T;
J2=1;
%弛豫步数
Thermal=1000;    
%采样数
sampling=1000;  
bins=sampling/20;

%采样间隔
interval=2;       

%预分配
m_sampling=zeros(1,sampling);
e_sampling=zeros(1,sampling);
cluster_max=zeros(1,sampling);
m=zeros(1,length(beta));
e=zeros(1,length(beta));
m2=zeros(1,length(beta));
cv=zeros(1,length(beta));
Variance=zeros(1,length(beta));
mean_cluster_max=zeros(1,length(beta));

for t=1:length(beta)
    %生成
    ising=generate_ising(n);
    %近邻
    nbor=neighbor(n,Nsite);
    for Relaxation_number=1:Thermal %弛豫
        [bonds]=configure_bonds(ising,nbor,n,J2,beta(t));
        [cluster]=search_cluster(bonds,nbor,Nsite,n);
        [cluster]=relabel(Nsite,cluster);
        %检查集团并画图
        %print_cluster(cluster,n,ising);
        [ising]=flip_spin(cluster,ising);
    end
    for sampling_number=1:sampling %采样
        for kk=1:interval
            [bonds]=configure_bonds(ising,nbor,n,J2,beta(t));
            [cluster]=search_cluster(bonds,nbor,Nsite,n);
            [cluster]=relabel(Nsite,cluster);
            [ising]=flip_spin(cluster,ising);
        end
        m_sampling(sampling_number)=abs(sum(sum(ising)))/Nsite;
        %[cluster_max(sampling_number)]=find_max_cluster(cluster);%找出最多的元素数量
        [e_sampling(sampling_number)]=calculate_energy(ising,n,nbor,J2);
    end
    
    mean_cluster_max(t)=mean(cluster_max);
    m(t)=mean(m_sampling);
    e(t)=mean(e_sampling);
    m2(t)=beta(t)*Nsite*std(m_sampling)^2;
    cv(t)=beta(t)^2 * n^2 * std(e_sampling)^2;
    Variance(t)=error_bar(e_sampling,bins);
    
    fprintf('temperature is %f\n',T(t));
end

% talbe=tabulate(energ(:))
% x=talbe(:,1); y=talbe(:,2);
% bar(x,y,0.4)-w 

figure(1);
subplot(2,2,1);errorbar(T,e,Variance,'k>:'); xlabel('T'); ylabel('E');  title('L=4,采样20000');  hold on
subplot(2,2,2);plot(T,m,'bo');   xlabel('T'); ylabel('M');     hold on
subplot(2,2,3);plot(T,cv,'bo');  xlabel('T'); ylabel('C_{v}'); hold on
subplot(2,2,4);plot(T,m2,'bo');  xlabel('T'); ylabel('\chi');  hold on

% data=[T',e',m',cv',m2',Variance'];
% save C:\Users\LENOVO\Desktop\e6.dat -ascii  data
% h=figure(1);
% saveas(h,'C:\Users\LENOVO\Desktop\figure_6.pdf');


%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
function print_cluster(cluster,n,ising)%画出cluster图
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

%--------------------------------------------------------------------------
function [cluster_max]=find_max_cluster(cluster)%找出最多的元素数量
cluster=cluster(:);
table = tabulate(cluster);
[~,idx] = max(table(:,2));%[maxCount,idx] = max(table(:,2));
%获取出现次数最多的元素
cluster_max=sum(cluster(:)==table(idx));

%--------------------------------------------------------------------------
%generate fe ising model
function [ising]=generate_ising(n)
ising=ones(1,n*n);
for i=1:n^2
    if rand>0.5
        ising(i)=-ising(i);
    end
end

%--------------------------------------------------------------------------
%generate table of near neighbor
function [nbor]=neighbor(n,Nsite)
nbor=zeros(Nsite,4);
for ispin=1:Nsite
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

%--------------------------------------------------------------------------
%configure bonds
function [bonds]=configure_bonds(ising,nbor,n,J2,beta)
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

%--------------------------------------------------------------------------
%apply the H-K algorithm and search cluster with PBC
function [cluster]=search_cluster(bonds,nbor,Nsite,n)
cluster=zeros(1,n*n);
largest_label=0;
label=(1:Nsite);
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
%find root label
function [y,label]=find_root(m,label)
y=m;
while label(y)~=y
    y=label(y);
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
%search nbor site
function sitenbor=searchnbor(number,nbor,cen)
if number==1
    sitenbor=nbor(cen,1);
elseif number==2
    sitenbor=nbor(cen,2);
% elseif number==3
%     sitenbor=nbor(nbor(cen,2),1);
% elseif number==4
%     sitenbor=nbor(nbor(cen,2),3);
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
%relabel cluster
function [cluster]=relabel(Nsite,cluster)
relabel_number=1;
for i=1:Nsite
    if cluster(cluster==i)~=0
        cluster(cluster==i)=relabel_number;
        relabel_number=relabel_number+1;
    end
end

%--------------------------------------------------------------------------
%flip cluster's spins
function [ising]=flip_spin(cluster,ising)
for i=1:max(max(cluster))
    if rand>0.5 %成功翻转
        [n]=find(cluster==i);
        for j=1:length(n)
            ising(n(j))=-ising(n(j));
        end
    end
end

%--------------------------------------------------------------------------
%calculate
function [energy]=calculate_energy(ising,n,nbor,J2)
energy=0.0;
for i=1:n*n    
    right=nbor(i,1);
    up=nbor(i,2);
    energy=energy - J2 * ising(i) * ( ising(up) + ising(right) );
end
energy =energy/n^2;
