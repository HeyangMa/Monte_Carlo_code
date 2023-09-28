% Parallel Tempering
% Ising model
% author MHY
% 2022-08-01
%% 参数
n = 4;
T = [.1, .33];
alpha = .5;
iter = 10000;%迭代数
sampler=5000;
for i=1:20
    [e_mean(i)]=calculate_one_temp(n,T,iter,sampler,alpha);
    fprintf('temperature is %f\t',T); fprintf('--已完成--\n');
    T=T+0.2;
end
figure(2);hold on
plot((0.1:0.2:0.1+19*0.2),e_mean,'r*-')
%% main
function [e_mean]=calculate_one_temp(n,T,iter,sampler,alpha)
nbor = neighbor(n);
ising = ones(length(T), n*n);
for a=1:iter
    [ising]=mcmcmc(alpha,ising,T,nbor);
end
for a=1:sampler
    for j=1:20
        [ising]=mcmcmc(alpha,ising,T,nbor);
    end
    e_samp(a)=Energy(ising(1,:),nbor);
end
e_mean=mean(e_samp);
end
%%
function [ising]=mcmcmc(alpha,ising,T,nbor)
if rand(1) < alpha
    % Normal Gibbs updating: parallel step
    [ising] = updateAllStates(ising,T,nbor);
else
    % Propose to exchange two states
    i = randsample(length(T)-1,1);
    DeltaE = Energy(ising(i,:),nbor)-Energy(ising(i+1,:),nbor);
    if rand(1) < min( 1,exp(DeltaE*( 1/T(i) - 1/T(i+1) )) )
        aux1 = ising(i,:);
        ising(i,:) = ising(i+1,:);
        ising(i+1,:) = aux1;
    end
end
end
%%
function [ising] = updateAllStates(ising,T,nbor)
for i=1:length(T)
    [ising(i,:)] = updateGibbs(ising(i,:),T(i),nbor);
end
end
%% 一个 monte carlo 步
function [ising] = updateGibbs(ising,T,nbor)
for i=1:length(ising)
    %rand pick one site from ising
    csite=unidrnd(length(ising));
    right=nbor(csite,1); up=nbor(csite,2);  left=nbor(csite,3);  down=nbor(csite,4);
    %calculate E_old
    E_old=-ising(csite)*(ising(up)+ising(right)+ising(down)+ising(left));
    %calculate E_new
    E_new=ising(csite)*(ising(up)+ising(right)+ising(down)+ising(left));
    % 能量差
    if rand<exp(-1*(E_new-E_old)/T)%exp(-1*E_new/T)/(exp(-1*E_old/T)+exp(-1*E_new/T))%
        ising(csite)=-ising(csite);
    end
end
end
%% 物理量
function [energy]=Energy(ising,nbor)
energy = 0.0;
for j = 1:length(ising)
    right=nbor(j,1); up=nbor(j,2);
    energy=energy - ising(j)*ising(right);           %right
    energy=energy - ising(j)*ising(up);              %up
end
energy = energy / length(ising);
end
%% 邻居
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