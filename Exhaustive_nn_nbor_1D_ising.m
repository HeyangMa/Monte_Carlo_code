%穷举: many-body interaction
%H=\sum -J2*SiSj
% close all
clear
%% 改变n,变为不同尺寸
n=16;
J2=1;
%数组大小
ising=-ones(1,n);
nbor=neighbor_pbc(n);
energy=zeros(2^(n),1);
m=zeros(2^(n),1);
%初始全为-1,计算
[energy(1),m(1)]=calculate(ising,nbor,n,J2);
%2^n-1个值
for i=1:length(energy)-1
    a=str2double(dec2bin(i));
    a=num2str(a,append('%0',num2str(n),'d'));
    for j=1:n
        if str2double(a(j)) == 1
            ising(j)=1;
        else
            ising(j)=-1;
        end
    end
    [energy(i+1),m(i+1)]=calculate(ising,nbor,n,J2);
end
%根据公式求不同温度下的能量期望
t=(0:0.1:4);
e_mean=zeros(1,length(t));
m_mean=zeros(1,length(t));
for i=1:length(t)
    %计算当前温度下的能量期望
    e_mean(i)=sum(energy .* exp(-(1/t(i))*n*energy)) / sum(exp(-(1/t(i))*n*energy));
    m_mean(i)=sum(m .* exp(-(1/t(i))*n*energy)) / sum(exp(-(1/t(i))*n*energy));
end
tabulate(energy)
figure(1);hold on;plot(t,e_mean,'k*-');xlabel('T');ylabel('E')
figure(2);hold on;plot(t,m_mean,'k*-');xlabel('T');ylabel('M')

%% generate neighbor
function [nbor]=neighbor_pbc(n)
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
%% calculate_energy
function [energy,m]=calculate(ising,nbor,n,J2)
energy=0.0;

for i=1:n
    
    energy=energy - J2 * ising(i) * ( ising(nbor(nbor(i,1),1))  );

end

energy=energy/(n);
m=abs(sum(ising))/n;
end