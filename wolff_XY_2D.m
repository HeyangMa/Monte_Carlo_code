%Author: Ma HY
%Date: 2023-04
% H =\sum -J * cos(\theta_1 \theta_2)
%% --------参数设置--------
n=6;               %维度格点数
J2=1;              %耦合常数J
T=(0:0.1:2);       %温度设置
Thermal=1000;      %弛豫步数
bins=100;           %bins数目
bsteps=1000;       %每个bin内的步数
%%
for i = 1:length(T)
    [e(i),m(i)]=mcmc(n, T(i), Thermal, bins, bsteps, J2);
end
figure(1)
plot(T,e,'ro')


%% a certain temperature Monte Carlo simulate
function [e_mean,m_mean]=mcmc(n, T, Thermal, bins, bsteps, J2)
beta = 1 / T;
xy = initspin(n);
nbor = neighbor(n);
for j = 1:Thermal
    xy = one_step_MonteCarlo(xy, nbor, n, beta, J2);
end
for j = 1:bins*bsteps
    xy = one_step_MonteCarlo(xy, nbor, n, beta, J2);
    [e(j),m(j)] = calculate(xy, nbor, n, J2);
end
e_mean=mean(e);
m_mean=mean(m);
fprintf('temperature is %f\t',T); fprintf('--已完成--\n');
end

%% spin initialization
function [spin]=initspin(n)
spin = zeros(n*n,1);
end
%% neighbor relationship
function [neigh]=neighbor(n)
neigh = zeros( n * n, 4);
for ispin = 1:n*n
    iy= fix((ispin-1)/n)+1;
    ix=ispin-(iy-1)*n;
    ixp=ix+1-fix(ix/n)*n;
    iyp=iy+1-fix(iy/n)*n;
    ixm=ix-1+fix((n-ix+1)/n)*n;
    iym=iy-1+fix((n-iy+1)/n)*n;
    neigh(ispin,1)=(iy-1)*n+ixp;%右邻居
    neigh(ispin,4)=(iym-1)*n+ix;%下邻居
    neigh(ispin,3)=(iy-1)*n+ixm;%左邻居
    neigh(ispin,2)=(iyp-1)*n+ix;%上邻居
end
end
%% calculate delta{E}
function [delta_E]=calculate_delta_E(xy, xy_new, csite, nbor, nb, J2)
% H =\sum -J2 * cos(\theta_i_1 \theta_i_2)
E_old = -J2 * cos(xy(csite) - xy_new(nbor(csite, nb)));
E_new = -J2 * cos(xy_new(csite) - xy_new(nbor(csite, nb)));
delta_E = E_new - E_old;
end
%% one step Monte Carlo simulate
function [xy_new]=one_step_MonteCarlo(xy, nbor, n, beta, J2)
xy_new = xy;
%choose a random vector
alpha = pi * rand;
%pick a site as starting point
csite = unidrnd(n*n);
cluster = zeros(n * n, 1);
flag = zeros(n * n, 1);
%表示把current_site放入了需要查找的队列
flag(1) = csite;
xy_new(csite) = 2 * alpha - xy_new(csite);
%修改值为1，表示该格点意见被选择过
cluster(csite) = 1;
%表示待查找的点的个数+1，如果为1就结束
counter = 2;
%开始查找
while counter > 1
    %从最后一个被加入的点开始查找
    counter = counter - 1;
    csite = flag(counter);
    flag(counter) = 0;
    for i = 1:4
        %查找每一个近邻
        next_site = nbor(csite, i);
        %需要满足条件
        deltaE = calculate_delta_E(xy, xy_new, csite, nbor, i, J2);
        %P=e^(-\beta δ{E})  therefore,
        probability = 1 - exp(-1 * beta * deltaE);
        if rand < probability && cluster(next_site) == 0
            xy_new(next_site) = 2 * alpha - xy_new(next_site);
            if xy_new(next_site) > 2 * pi
                xy_new(next_site) = xy_new(next_site) - 2 * pi;
            elseif xy_new(next_site) < 0
                xy_new(next_site) = xy_new(next_site) + 2 * pi;
            end
            %将满足条件的点放入待查找的队列中
            flag(counter) = next_site;
            cluster(next_site) = 1;
            counter = counter + 1;
        end
    end
end
end
%% calculate physical quantities in equilibrium
function [energy,m]=calculate(xy, nbor, n, J2)
energy = 0.0;
for j = 1:n*n
    up = nbor(j, 2);
    left = nbor(j, 3);
    energy = energy - J2 * (cos(xy(j) - xy(left)));
    energy = energy - J2 * (cos(xy(j) - xy(up)));
end
energy = energy / (n * n);
m=sqrt(sum(cos(xy))^2+sum(sin(xy))^2)/ (n * n);
end