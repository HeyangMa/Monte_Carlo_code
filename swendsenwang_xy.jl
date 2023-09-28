"""
示例的主文件 IXYSA_wolff.jl
Author: Ma HY
Date: 2023-05
"""
# H=-\sum J*cos(\theta_{i} - \theta_{j})

#xy initialization
function init_xy(n::Int)
    spin = zeros(Float64, n * n, 1)
    return spin
end

#neighbor relationship
function neighbor(ll)
    neigh = Array{Int}(undef, ll * ll, 4)
    for s0 = 1:ll*ll
        x0 = mod(s0 - 1, ll)
        y0 = div(s0 - 1, ll)
        x1 = mod(x0 + 1, ll)
        x2 = mod(x0 - 1 + ll, ll)
        y1 = mod(y0 + 1, ll)
        y2 = mod(y0 - 1 + ll, ll)
        neigh[s0, 1] = 1 + x1 + y0 * ll
        neigh[s0, 2] = 1 + x0 + y1 * ll
        neigh[s0, 3] = 1 + x2 + y0 * ll
        neigh[s0, 4] = 1 + x0 + y2 * ll
    end
    return neigh

end

#[c,d]=project_t(a,b)
function project_t(a, b)
    theta = a
    theta_proj = b
    # 计算角向量和投影方向向量
    vec_angle = [cos(theta), sin(theta)]
    proj_dir = [cos(theta_proj), sin(theta_proj)]
    # 计算投影值
    proj_value = vec_angle[1] * proj_dir[1] + vec_angle[2] * proj_dir[2]

    ## --------------给定两个角度和方向---------------------
    theta1 = a
    theta2 = b

    # 计算两个方向向量
    dir1 = [cos(theta1), sin(theta1)]
    dir2 = [cos(theta2), sin(theta2)]

    # 计算反射向量
    reflection_dir = dir1 - (2 * (dir1[1] * dir2[1] + dir1[2] * dir2[2])) * dir2

    # 计算反射角
    theta_refl = atan(reflection_dir[2], reflection_dir[1])
    if theta_refl < 0
        theta_refl = 2 * pi + theta_refl # 保证反射角度在0到360之间
    end
    angle_refl = theta_refl * 180 / pi
    d = angle_refl * pi / 180

    proj_value, d
end

#configure bonds
function configure_bonds(xy, nbor, n, beta, J, alpha)
    bonds = zeros(Int, n * n, 2)
    #计算投影
    xy_pro = zeros(Float64, n * n, 1)
    for j = 1:n*n
        xy_pro[j], ~ = project_t(xy[j], alpha)
    end
    #place bonds
    for j = 1:n*n
        if rand() < max(0, 1 - exp(-2 * J * beta * xy_pro[j] * xy_pro[nbor[j, 1]]))
            bonds[j, 1] = 1
        end
        if rand() < max(0, 1 - exp(-2 * J * beta * xy_pro[j] * xy_pro[nbor[j, 2]]))
            bonds[j, 2] = 1
        end
    end
    #----End----
    return bonds
end

#search nbor site
function searchnbor(number, nbor, cen)
    if number == 1
        sitenbor = nbor[cen, 1]
    elseif number == 2
        sitenbor = nbor[cen, 2]
        # elseif number == 3
        #     sitenbor = nbor[nbor[cen, 2], 1]
        # elseif number == 4
        #     sitenbor = nbor[nbor[cen, 2], 3]
        # elseif number==5
        #     sitenbor=nbor(nbor(nbor(cen,3),3),2);
        # elseif number==6
        #     sitenbor=nbor(nbor(nbor(cen,3),2),2);
        # elseif number==7
        #     sitenbor=nbor(nbor(nbor(cen,1),2),2);
        # elseif number==8
        #     sitenbor=nbor(nbor(nbor(cen,1),1),2);
        # elseif number==9
        #     sitenbor=nbor(nbor(cen,2),2);
        # elseif number==10
        #     sitenbor=nbor(nbor(cen,3),3);
    end
    return sitenbor
end

#find root label
function find_root(m, label)
    y = m
    while label[y] != y
        y = label[y]
    end
    y, label
end

#union bonds
function union(notzerolabellist, label)
    var = Array{Int}(undef, length(notzerolabellist), 1)
    for k = 1:length(notzerolabellist)
        aaa = find_root(notzerolabellist[k], label)
        var[k] = aaa[1]
        label = aaa[2]
    end
    a::Int64 = findmin(var)[1]
    for k = 1:length(notzerolabellist)
        label[var[k]] = a
    end
    m = a
    label, m
end

#apply the H-K algorithm and search cluster with PBC
function search_cluster(bonds, nbor, n)
    cluster = zeros(Int64, n * n, 1) #cluster label
    label = Array{Int64}(1:n*n) #root label
    largest_label = 0
    for j = 1:n*n
        #1:right  2:up 
        bondlist = [bonds[j, 1], bonds[j, 2]]
        bondsum = sum(bondlist)
        if cluster[j] == 0 #current site unlabeled
            if bondsum == 0 #no bond connection
                largest_label = largest_label + 1
                cluster[j] = largest_label
            elseif bondsum == 1 ##one bond connection
                number = findall(==(1), bondlist)
                sitenbor = searchnbor(number[1], nbor, j)
                if cluster[sitenbor] == 0
                    largest_label = largest_label + 1
                    cluster[j] = largest_label
                    cluster[sitenbor] = largest_label
                else
                    aaa = find_root(cluster[sitenbor], label)
                    cluster[j] = aaa[1]
                    label = aaa[2]
                end
            elseif bondsum > 1 #more than one bond connection
                numberlist = findall(==(1), bondlist)
                #use the following array to memory
                varlist = zeros(Int, length(numberlist), 1)   #memory cluster label
                sitenbor1 = zeros(Int, length(numberlist), 1) #memory site  label
                for k = 1:length(numberlist)
                    sitenbor1[k] = searchnbor(numberlist[k], nbor, j)
                    varlist[k] = cluster[sitenbor1[k]]
                end
                labellist = varlist[findall(!=(0), varlist)]
                if labellist == []
                    largest_label = largest_label + 1
                    for k = 1:length(numberlist)
                        cluster[sitenbor1[k]] = largest_label
                    end
                    cluster[j] = largest_label
                else
                    aaa = union(labellist, label)
                    label = aaa[1]
                    minlabel = aaa[2]
                    for k = 1:length(numberlist)
                        cluster[sitenbor1[k]] = minlabel
                    end
                    cluster[j] = minlabel
                end
            end
        else #The current site has a non-zeros label
            if bondsum == 0 #no bond
                continue
            elseif bondsum > 0 #more than one bond (current site bring one)
                numberlist = findall(==(1), bondlist)
                #use the following array to memory
                varlist = zeros(Int, length(numberlist), 1)   #memory cluster label
                sitenbor1 = zeros(Int, length(numberlist), 1) #memory site label
                for k = 1:length(numberlist)
                    sitenbor1[k] = searchnbor(numberlist[k], nbor, j)
                    varlist[k] = cluster[sitenbor1[k]]
                end
                labellist = varlist[findall(!=(0), varlist)]
                if labellist == []
                    for k = 1:length(numberlist)
                        aaa = find_root(cluster[j], label)
                        a = aaa[1]
                        label = aaa[2]
                        cluster[sitenbor1[k]] = a
                    end
                else
                    aaa = union(labellist, label)
                    label = aaa[1]
                    minlabel = aaa[2]
                    aaa = find_root(cluster[j], label)
                    a = aaa[1]
                    label = aaa[2]
                    sminlabel = findmin([minlabel, a])[1]
                    label[minlabel] = sminlabel
                    label[a] = sminlabel
                    for k = 1:length(numberlist)
                        cluster[sitenbor1[k]] = minlabel
                    end
                    cluster[j] = minlabel
                end
            end
        end
    end
    for j = 1:n*n
        aaa = find_root(cluster[j], label)
        cluster[j] = aaa[1]
        label = aaa[2]
    end
    return cluster
end

#flip cluster's spins (this is one Monte Carlo step)
function flip_spin(cluster, xy, alpha)
    for i = 1:findmax(cluster)[1]
        if rand() > 0.5 && findall(==(i), cluster) != []  #成功翻转 && 有这一标签的集团
            n = findall(==(i), cluster)

            for j = 1:length(n)
                #xy
                ~, xy[n[j]] = project_t(xy[n[j]], alpha)
            end
        end
    end
    return xy
end

#calculate other phys quantities
function Jackknife(data)
    ave_data = zeros(Float64, length(data), 1)
    for i = 1:length(data)
        a = 0
        for j = 1:length(data)
            if j != i
                a = a + data[j]
            end
        end
        ave_data[i] = a / (length(data) - 1)
    end
    expect = sum(data) / length(data) - (length(data) - 1) * (sum(ave_data) / length(data) - sum(data) / length(data))
    b = ave_data - ones(Float64, length(data), 1) * (sum(ave_data) / length(data))
    delta_expect = sqrt(((length(data) - 1) / length(data)) * sum(b .^ 2))
    expect, delta_expect
end

#sample variance
function var(data)
    mean_data = sum(data) / length(data)
    s = 0
    for i = 1:length(data)
        s = s + (data[i] - mean_data)^2
    end
    s = s / (length(data) - 1)
    return s
end

#count integer vortex
function count_vortex(xy, nbor, n)
    theta = zeros(Float64, 4, 1)
    v = zeros(Float64, 4, 1)
    vortex = 0
    error = 0.01
    for i = 1:n^2
        theta[1] = xy[nbor[nbor[i, 1], 2]]
        theta[2] = xy[nbor[i, 2]]
        theta[3] = xy[i]
        theta[4] = xy[nbor[i, 1]]
        v[1] = theta[2] - theta[1]
        v[2] = theta[3] - theta[2]
        v[3] = theta[4] - theta[3]
        v[4] = theta[1] - theta[4]
        for k = 1:4
            if v[k] <= -pi
                v[k] = v[k] + 2 * pi
            end
            if v[k] >= pi
                v[k] = v[k] - 2 * pi
            end
        end
        kv = sum(v) / 2 / pi
        if floor(kv) == 1
            vortex = vortex + 1
        end
        if floor(kv) == -1
            vortex = vortex + 1
        end
    end
    return vortex
end

#process_bindata
function process_bin(a, n, T, bsteps)
    #energy, m_xy, m_xy2, theta, theta2, H_x, H_y, I_x, I_y,m_xy_x,m_xy_y
    #average
    e = sum(a[:, 1] / n^2) / length(a[:, 1])
    m_xy = sum(a[:, 2] / n^2) / length(a[:, 2])
    m_xy2 = sum(a[:, 3] / n^2) / length(a[:, 3])
    theta = sum(a[:, 4] / n^2) / length(a[:, 4])
    theta2 = sum(a[:, 5] / n^2) / length(a[:, 5])

    #calculate high-order derivative
    #----cv----
    cv = var(a[:, 1]) / T^2

    #----some ms----
    ms_xy = var(a[:, 2]) / T
    ms_xy2 = var(a[:, 3]) / T

    #----fluctuation vortex----
    ms_theta = var(a[:, 4])
    ms_theta2 = var(a[:, 5])

    #----some Binder cumulation----
    R2_xy = (sum((a[:, 2] / n^2) .^ 4) / bsteps) / (sum((a[:, 2] / n^2) .^ 2) / bsteps)^2
    Br_xy = 2 * (1 - R2_xy / 2)
    #
    R2_xy2 = (sum((a[:, 3] / n^2) .^ 4) / bsteps) / (sum((a[:, 3] / n^2) .^ 2) / bsteps)^2
    Br_xy2 = 2 * (1 - R2_xy2 / 2)

    #rho
    H_x = sum(a[:, 6]) / length(a[:, 6])
    H_y = sum(a[:, 7]) / length(a[:, 7])
    I_x = sum(a[:, 8] .^ 2) / length(a[:, 8])
    I_y = sum(a[:, 9] .^ 2) / length(a[:, 9])
    rho_x = (H_x - I_x / T) / n^2
    rho_y = (H_y - I_y / T) / n^2
    rho = (rho_x + rho_y) / 2

    #chi_3
    variable1 = (sum((a[:, 10] / n^2) .^ 4) / length(a[:, 10])) - 3 * (sum((a[:, 10] / n^2) .^ 2) / length(a[:, 10]))^2
    variable2 = (sum((a[:, 11] / n^2) .^ 4) / length(a[:, 11])) - 3 * (sum((a[:, 11] / n^2) .^ 2) / length(a[:, 11]))^2
    chi_3 = (variable1 + variable2) * n^2 / T^3 / 3

    e, m_xy, m_xy2, theta, theta2, cv, ms_xy, ms_xy2, ms_theta, ms_theta2, Br_xy, Br_xy2, rho, chi_3
end

#calculate physical quantities in equilibrium
function calculate(xy, nbor, n, J)
    energy = 0.0
    for j = 1:n*n
        energy = energy - J * (cos(xy[j] - xy[nbor[j, 1]]))
        energy = energy - J * (cos(xy[j] - xy[nbor[j, 2]]))
    end

    # xyxy
    xy2 = copy(xy) * 2
    for i = 1:n^2
        while xy2[i] > 2 * pi || xy2[i] < 0
            if xy2[i] > 2 * pi
                xy2[i] = xy2[i] - 2 * pi
            end
            if xy2[i] < 0
                xy2[i] = xy2[i] + 2 * pi
            end
        end
    end
    m_xy = sqrt((sum(cos.(xy)))^2 + (sum(sin.(xy)))^2)
    m_xy2 = sqrt((sum(cos.(xy2)))^2 + (sum(sin.(xy2)))^2)
    theta = count_vortex(xy, nbor, n)
    theta2 = count_vortex(xy2, nbor, n)

    #I_x I_y
    H_x = 0
    H_y = 0
    I_x = 0
    I_y = 0
    for j = 1:n^2
        H_x = H_x + J * (cos(xy[j] - xy[nbor[j, 1]]))
        H_y = H_y + J * (cos(xy[j] - xy[nbor[j, 2]]))
        I_x = I_x + J * (sin(xy[j] - xy[nbor[j, 1]]))
        I_y = I_y + J * (sin(xy[j] - xy[nbor[j, 2]]))
    end

    #m_xy_x m_xy_y
    m_xy_x = sum(cos.(xy)) / n^2
    m_xy_y = sum(sin.(xy)) / n^2
    energy, m_xy, m_xy2, theta, theta2, H_x, H_y, I_x, I_y, m_xy_x, m_xy_y
end

#a certain temperature Monte Carlo simulate
function mcmc(n, T, Thermal, bins, bsteps, J, certain, diff)
    beta = 1 / T
    xy = init_xy(n)
    nbor = neighbor(n)
    for j = 1:Thermal
        # one Monte Carlo step
        alpha = 2 * pi * rand()
        bonds = configure_bonds(xy, nbor, n, beta, J, alpha)
        cluster = search_cluster(bonds, nbor, n)
        xy = flip_spin(cluster, xy, alpha)
    end
    # certain T
    if certain == 1
        write_data(n, T, beta, xy, nbor, bins, bsteps, J)
    end
    # different T
    if diff == 1
        write_data2(n, T, beta, xy, nbor, bins, bsteps, J)
    end
end

#write data at certain T
function write_data(n, T, beta, xy, nbor, bins, bsteps, J)
    mdata = zeros(Float64, bins * bsteps, 5)
    dizhi = string(k, "Tc_wolff_IXYSA_L", n, ".txt")
    f = open(dizhi, "a")
    for j = 1:bins*bsteps
        # one Monte Carlo step
        alpha = 2 * pi * rand()
        bonds = configure_bonds(xy, nbor, n, beta, J, alpha)
        cluster = search_cluster(bonds, nbor, n)
        xy = flip_spin(cluster, xy, alpha)
        mdata[j, 1], mdata[j, 2], mdata[j, 3], mdata[j, 4], mdata[j, 5] = calculate(ising, xy, nbor, n, A, B, C)
        println(f, mdata[i, 1], "    ", mdata[i, 2], "    ", mdata[i, 3])
    end
    close(f)
end

#write data at different T with error
function write_data2(n, T, beta, xy, nbor, bins, bsteps, J)
    bin = zeros(Float64, bins, 14)
    for j = 1:bins
        bstep = zeros(Float64, bsteps, 11)
        for i = 1:bsteps
            # one Monte Carlo step
            alpha = 2 * pi * rand()
            bonds = configure_bonds(xy, nbor, n, beta, J, alpha)
            cluster = search_cluster(bonds, nbor, n)
            xy = flip_spin(cluster, xy, alpha)
            # measure
            #energy,     m_xy,        m_xy2,       theta,       theta2,      H_x,         H_y,         I_x,         I_y      m_xy_x     m_xy_y
            bstep[i, 1], bstep[i, 2], bstep[i, 3], bstep[i, 4], bstep[i, 5], bstep[i, 6], bstep[i, 7], bstep[i, 8], bstep[i, 9], bstep[i, 10], bstep[i, 11] = calculate(xy, nbor, n, J)
        end
        #process bindata
        #e, m_xy, m_xy2, theta, theta2, cv, ms_xy, ms_xy2, ms_theta, ms_theta2, Br_xy, Br_xy2, rho, chi_3
        bin[j, 1], bin[j, 2], bin[j, 3], bin[j, 4], bin[j, 5], bin[j, 6], bin[j, 7], bin[j, 8], bin[j, 9], bin[j, 10], bin[j, 11], bin[j, 12], bin[j, 13], bin[j, 14] = process_bin(bstep, n, T, bsteps)
    end

    #calculate average with error bar
    e = sum(bin[:, 1]) / bins
    e_error = sqrt(var(bin[:, 1]))
    #
    m_xy = sum(bin[:, 2]) / bins
    m_xy_error = sqrt(var(bin[:, 2]))
    #
    m_xy2 = sum(bin[:, 3]) / bins
    m_xy2_error = sqrt(var(bin[:, 3]))
    #
    theta = sum(bin[:, 4]) / bins
    theta_error = sqrt(var(bin[:, 4]))
    #
    theta2 = sum(bin[:, 5]) / bins
    theta2_error = sqrt(var(bin[:, 5]))

    #Jackknife: cv, ms_xy, ms_xy2, ms_theta, ms_theta2, Br_xy, Br_xy2, rho, chi_3
    cv, cv_error = Jackknife(bin[:, 6])
    ms_xy, ms_xy_error = Jackknife(bin[:, 7])
    ms_xy2, ms_xy2_error = Jackknife(bin[:, 8])
    ms_theta, ms_theta_error = Jackknife(bin[:, 9])
    ms_theta2, ms_theta2_error = Jackknife(bin[:, 10])
    Br_xy, Br_xy_error = Jackknife(bin[:, 11])
    Br_xy2, Br_xy2_error = Jackknife(bin[:, 12])
    rho, rho_error = Jackknife(bin[:, 13])
    chi_3, chi_3_error = Jackknife(bin[:, 14])

    #print data 
    dizhi = string("C://Users//MHY//Desktop//wolff_XY_phyL", n, ".txt")
    f = open(dizhi, "a")
    println(f, T, "    ", e, "    ", e_error, "    ", m_xy, "    ", m_xy_error, "    ", m_xy2, "    ", m_xy2_error, "    ", theta, "    ", theta_error, "    ", theta2, "    ", theta2_error, "    ", cv, "   ", cv_error, "    ", ms_xy, "    ", ms_xy_error, "    ", ms_xy2, "    ", ms_xy2_error, "    ", ms_theta, "    ", ms_theta_error, "    ", ms_theta2, "    ", ms_theta2_error, "    ", Br_xy, "    ", Br_xy_error, "    ", Br_xy2, "    ", Br_xy2_error, "    ", rho, "    ", rho_error, "    ", chi_3, "    ", chi_3_error)
    close(f)
end

#--------参数设置-------- 
# txt1 = readline("read_n.in")
# n = parse(Int64, txt1)
n = [10]
J = 1
T = (0:0.1:2)
Thermal = 5000       #弛豫步数
bins = 100            #bins数目
bsteps = 100        #每个bin内的步数

#--------控制--------
certain = 0
diff = 1
#--------------------------------------------------------------------------
@time begin
    for j = 1:length(n)
        for i = 1:length(T)
            mcmc(Int(n[j]), T[i], Thermal, bins, bsteps, J, certain, diff)
        end
    end
end