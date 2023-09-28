%metropolis: no_interaction_ising
%H= - \sum J*SiSj
%--------------------------------------------------------------------------
function metropolis_Ising_2D_no_interaction
n=10;
Nsite=n*n;
T=(0:0.2:5);
J=-1;
Thermal=1000;         %弛豫步数
sampling=1000;        %采样数
interval=100;         %采样间隔
%预分配
m_sampling=zeros(1,sampling);
e_sampling=zeros(1,sampling);
m=zeros(1,length(T));
e=zeros(1,length(T));
cv=zeros(1,length(T));
for t=1:length(T)
    
    beta=1.0/T(t);
    
    [site_ising]=generate_ising(n);
    
    for j=1:Thermal%弛豫
        [site_ising]=flip_spin(site_ising,n,Nsite,J,beta);
    end
    
    for j=1:sampling%采样
        for jj=1:interval
            [site_ising]=flip_spin(site_ising,n,Nsite,J,beta);
        end
        m_sampling(j)=abs(mean(mean(site_ising)));
        [e_sampling(j)]=calculate_energy(site_ising,n,J);
    end
    %average
    m(t)=mean(m_sampling);
    e(t)=mean(e_sampling);
    cv(t)=beta^2 * n^2 * std(e_sampling)^2;
    
    fprintf('temperature is %f\n',T(t));
end

figure(1);
plot(T,e,'b*')
xlabel('T')
ylabel('E')
hold on
figure(2);
plot(T,m,'b*')
xlabel('T')
ylabel('M')
hold on
figure(3);
plot(T,cv,'b*')
xlabel('T')
ylabel('C_{v}')
hold on

%generate ising model
function [site_ising]=generate_ising(n)
site_ising=ones(n,n);

%flip spins
function [site_ising]=flip_spin(site_ising,n,Nsite,J,beta)
%当前
unid_site=unidrnd(Nsite);
iy0=fix((unid_site-1)/n)+1;
ix0=unid_site-(iy0-1)*n;
%计算能量差
E_old=-1 * J * site_ising(iy0,ix0);
E_new=1 * J * site_ising(iy0,ix0);
deE=E_new-E_old;
if rand<exp(-deE*beta)
    site_ising(iy0,ix0)=-site_ising(iy0,ix0);
end

%calculate_energy
function [e_sampling]=calculate_energy(site_ising,n,J)
energy=0.0;
for j=1:n
    for i=1:n
        energy=energy -1 * J * site_ising(j,i);
    end
end
e_sampling=energy/n^2;