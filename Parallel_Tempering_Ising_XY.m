% Parallel Tempering
% xy-ising model
% author MHY
% 2022-08-01
%% 参数
n = 4;
T = [.1, .33, .40];
alpha = .5;
iter = 100000;%迭代数
sampler=5000;
for i=1:20
[e_mean(i)]=calculate_one_tem(n,T,iter,sampler,alpha);
T=T+0.2;
end
figure(2);hold on
plot((0.1:0.2:0.1+19*0.2),e_mean,'ro-')
%% main
function [e_mean]=calculate_one_tem(n,T,iter,sampler,alpha)
nbor = neighbor(n);
ising = ones(length(T), n*n);
xy = zeros(length(T), n*n);

for a=1:iter
    if rand(1) < alpha
        % Normal Gibbs updating: parallel step
        [ising,xy] = updateAllStates(ising,xy,T,nbor);
    else
        % Propose to exchange two states
        i = randsample(length(T)-1,1);
        DeltaE = Energy(ising(i+1,:),xy(i+1,:),nbor)-Energy(ising(i,:),xy(i,:),nbor);
        if rand(1) < min(1,exp(DeltaE/(1/T(i+1) - 1/T(i))))
            aux1 = ising(i,:);
            ising(i,:) = ising(i+1,:);
            ising(i+1,:) = aux1;
            aux2 = xy(i,:);
            xy(i,:) = xy(i+1,:);
            xy(i+1,:) = aux2;
        end
    end
end
for a=1:sampler
    if rand(1) < alpha
        % Normal Gibbs updating: parallel step
        [ising,xy] = updateAllStates(ising,xy,T,nbor);
    else
        % Propose to exchange two states
        i = randsample(length(T)-1,1);
        DeltaE = Energy(ising(i+1,:),xy(i+1,:),nbor)-Energy(ising(i,:),xy(i,:),nbor);
        if rand(1) < min(1,exp(DeltaE/(1/T(i+1) - 1/T(i))))
            aux1 = ising(i,:);
            ising(i,:) = ising(i+1,:);
            ising(i+1,:) = aux1;
            aux2 = xy(i,:);
            xy(i,:) = xy(i+1,:);
            xy(i+1,:) = aux2;
        end
    end
    e_samp(a)=Energy(ising(1,:),xy(1,:),nbor);
end
    e_mean=mean(e_samp);
end
%%
function [ising,xy] = updateAllStates(ising,xy, T,nbor)
for i=1:length(T)
    [ising(i,:),xy(i,:)] = updateGibbs(ising(i,:),xy(i,:), T(i),nbor);
end
end
%% 一个 monte carlo 步
function [ising,xy] = updateGibbs(ising,xy,T,nbor)
        if rand<=0.5
            %rand pick one site from xy
            csite=unidrnd(length(ising));
            right=nbor(csite,1); up=nbor(csite,2);  left=nbor(csite,3);  down=nbor(csite,4);
            %calculate energy (old)
            E_old=-cos(xy(csite)-xy(left));                                      %left
            E_old=E_old - cos(xy(csite)-xy(right));                              %right
            E_old=E_old - ising(csite)*ising(up)*cos(xy(csite)-xy(up));          %up
            E_old=E_old - ising(csite)*ising(down)*cos(xy(csite)-xy(down));      %down
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
            E_new=-cos(angle-xy(left));                                      %left
            E_new=E_new - cos(angle-xy(right));                              %right
            E_new=E_new - ising(csite)*ising(up)*cos(angle-xy(up));          %up
            E_new=E_new - ising(csite)*ising(down)*cos(angle-xy(down));      %down
            %calculate delta energy
            deltaE=E_new-E_old;
            if rand<exp(-1*deltaE/T)
                xy(csite)=angle;
            end
        else
            %rand pick one site from ising
            csite=unidrnd(length(ising));
            right=nbor(csite,1); up=nbor(csite,2);  left=nbor(csite,3);  down=nbor(csite,4);
            %calculate energy (old)
            E_old=-cos(xy(csite)-xy(left));                                      %left
            E_old=E_old - cos(xy(csite)-xy(right));                              %right
            E_old=E_old -ising(csite)*ising(up)*cos(xy(csite)-xy(up));          %up
            E_old=E_old - ising(csite)*ising(down)*cos(xy(csite)-xy(down));      %down
            %calculate energy (new)
            E_new=-cos(xy(csite)-xy(left));                                      %left
            E_new=E_new - cos(xy(csite)-xy(right));                              %right
            E_new=E_new + ising(csite)*ising(up)*cos(xy(csite)-xy(up));          %up
            E_new=E_new + ising(csite)*ising(down)*cos(xy(csite)-xy(down));      %down
            % 能量差
            deltaE=E_new-E_old;
            if rand<exp(-1*deltaE/T)
                ising(csite)=-ising(csite);
            end
        end
end
%% 计算物理量
function [energy]=Energy(ising,xy, nbor)
energy = 0.0;
for j = 1:length(ising)
    right=nbor(j,1); up=nbor(j,2);
    energy=energy - cos(xy(j)-xy(right));                              %right
    energy=energy - ising(j)*ising(up)*cos(xy(j)-xy(up));              %up
end
energy = energy / length(ising);
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