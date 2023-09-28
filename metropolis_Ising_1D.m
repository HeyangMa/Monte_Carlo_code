%metropolis: many-body interaction
%H=\sum -J2*SiSj 
%--------------------------------------------------------------------------
%parameter setting
function metropolis_Ising_1D

%晶格长度
n=10;

%温度
T=(0:0.1:4); beta=1./T;

%耦合常数
J2=1;

%弛豫步数
Thermal=5000;   
%采样数
sampling=5000;  
%采样间隔
interval=100;    

%预分配
m_sampling=zeros(1,sampling);
e_sampling=zeros(1,sampling);
m=zeros(1,length(T));
e=zeros(1,length(T));


%邻居结构
[nbor]=neighbor(n);

for t=1:1:length(beta)
    fprintf('temperature is %f\n',T(t));
    [ising]=generate(n);%生成
    for j=1:Thermal%弛豫
        [ising]=flip_spin(ising,nbor,n,beta(t),J2);
    end
    for j=1:sampling%采样
        for kk=1:interval
            [ising]=flip_spin(ising,nbor,n,beta(t),J2); 
        end
        %sample
        m_sampling(j)=abs(sum(ising))/n;
        [e_sampling(j)]=calculate(ising,nbor,n,J2);
    end
    %average
    e(t)=mean(e_sampling);
    m(t)=mean(m_sampling);
end

figure(1);
plot(T,e,'ro-'); xlabel('T'); ylabel('E');   hold on

figure(2);
plot(T,m,'ro-');   xlabel('T'); ylabel('M'); hold on



%--------------------------------------------------------------------------
function [ising]=generate(n)
ising=ones(n,1);

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

function [ising]=flip_spin(ising,nbor,n,beta,J2)

cen=unidrnd(n);%当前

deE=2 * J2 * ising(cen) * ( ising(nbor(cen,1)) + ising(nbor(cen,2)) ) ;

if rand<exp(-deE*beta)
    ising(cen)=-ising(cen);
end

function [energy]=calculate(ising,nbor,n,J2)
energy=0.0;

for i=1:n
    
    energy=energy - J2 * ising(i) * ( ising(nbor(i,1))  );

end

energy=energy/(n);