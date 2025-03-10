function m0 = ustc_litho(target_lat,target_lon)

% ak135 model
[z0, rho0, vp0, vs0, qk0, qm0] = ak135( 'cont' );

% USTClitho2.0 速度模型
% 读取 USTClitho2.0.txt 数据
data = load('./velocity_model/USTClitho2.0.txt');
vp = data(:, 4);  % P 波速度 Vp (km/s)
vs = data(:, 5);  % S 波速度 Vs (km/s)
depth = data(:, 3);
lon = data(:, 1);
lat = data(:, 2);

% 获取所有唯一的深度值
depths = unique(depth);
num_depths = length(depths);

% 预存储每个深度层的插值函数
Vpinp = cell(num_depths, 1);
Vsinp = Vpinp;
% 对每个深度层进行处理
for k = 1:num_depths
    current_depth = depths(k);
    % 提取当前深度的数据点
    mask = (depth == current_depth);
    current_lon = lon(mask);
    current_lat = lat(mask);
    current_vp = vp(mask);
    current_vs = vs(mask);
    
    Fvp = scatteredInterpolant(current_lon, current_lat, current_vp, 'linear', 'none');

    Fvs = scatteredInterpolant(current_lon, current_lat, current_vs, 'linear', 'none');

    Vpinp{k} = Fvp;
    Vsinp{k} = Fvs;

end

m0 = zeros(length(depths),4);
for i = 1:length(depths)
    Fvp = Vpinp{i};
    Fvs = Vsinp{i};

    vp_interp = Fvp(target_lon, target_lat);
    vs_interp = Fvs(target_lon, target_lat);

    m0(i,1) = depths(i);
    m0(i,2) = vp_interp;
    m0(i,3) = vs_interp;

end
m0(:,4) = 2.70;


%%-----------AK135 model----------------------
keepz= z0 >= max(depths);
m_pad = [z0(keepz) vp0(keepz) vs0(keepz) rho0(keepz)];
ztemp_raw = m_pad(:,1);
vp_raw = m_pad(:,2);
vs_raw = m_pad(:,3);
rho_raw = m_pad(:,4);
% 去重并保持顺序
[ztemp_unique, idx_unique] = unique(ztemp_raw, 'stable');
vp_unique = vp_raw(idx_unique);
vs_unique = vs_raw(idx_unique);
rho_unique = rho_raw(idx_unique);

m_pad = [ztemp_unique,vp_unique,vs_unique,rho_unique];

m0 = [m0;m_pad];

return;
