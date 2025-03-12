function gridStruct = getVelocityModel(ModelType, gridStruct, npts)
%GETVELOCITYMODEL Generate velocity model based on specified type and parameters
%
% Inputs:
%   ModelType  - String, type of velocity model: '1D', '2D', or '3D'
%   gridStruct - Structure containing grid information:
%                For 1D: nx, z
%                For 3D: x, y, z, LatMin, LatMax, LonMin, LonMax
%   npts       - Number of points for interpolation (used in 3D model)
%
% Outputs:
%   gridStruct - Updated structure with velocity information:
%                For 1D: adds vp, vs fields
%                For 3D: adds VP, VS, Fvp, Fvs fields
%
% Note: 2D implementation is currently placeholder only

% Input validation
if ~ischar(ModelType) || ~ismember(ModelType, {'1D', '2D', '3D'})
    error('ModelType must be one of: 1D, 2D, or 3D');
end

gridStruct.ModelType = ModelType;
switch ModelType
    case '1D'
        % Get 1D velocity model from CRUST1.0
        model = obtain_crust1_QB();
        nx = gridStruct.nx;
        
        % Interpolate P and S wave velocities to grid points
        vp = interp1(model(:,1), model(:,2), gridStruct.z, 'nearest', 'extrap');
        vs = interp1(model(:,1), model(:,3), gridStruct.z, 'nearest', 'extrap');
        vel = repmat(vp(:), 1, nx);
        vel_s = repmat(vs(:), 1, nx);
        
        % Apply smoothing to velocity models using moving average
        N = 5;  % Smoothing window size
        [vel,~] = moving_avg(vel, N, 'constant', 2);
        [vel,~] = moving_avg(vel, N, 'constant');
        [vel_s,~] = moving_avg(vel_s, N, 'constant', 2);
        [vel_s,~] = moving_avg(vel_s, N, 'constant');
        
        gridStruct.vp = vel;
        gridStruct.vs = vel_s;
        
    case '2D'
        % TODO: Implement 2D velocity model
        warning('2D velocity model not implemented yet');
        
    case '3D'
        % Extract geographical boundaries
        LatMin = gridStruct.LatMin;
        LatMax = gridStruct.LatMax;
        LonMin = gridStruct.LonMin;
        LonMax = gridStruct.LonMax;

        % Get AK135 reference model
        [z, rho, vp, vs, qk, qm] = ak135('cont');
        zmax = 800;
        dz = 0.5;

        % Initialize arrays for storing 3D model data
        X = []; Y = []; Z = []; VP = []; VS = [];
        knode = 0;
        
        % Create lat/lon grid for interpolation
        latall = linspace(LatMin, LatMax, npts);
        lonall = linspace(LonMin, LonMax, npts);
        
        % Loop through each grid point to build 3D model
        for i = 1:length(latall)
            for j = 1:length(lonall)
                knode = knode + 1;
                disp(['Processing node ', num2str(knode), ' of ', num2str(npts^2)]);
                
                % Get velocity model at current location
                lat = latall(i);
                lon = lonall(j);
                % note that m0 starts at 0 km depth. The elevation infomration is
                % missing, to include topography set if_topo to 1
                % if_topo = 1;
                % [m0,~] = obtain_crust1_v2(lat,lon,[],if_topo);
                m0 = obtain_ustc_litho(lat,lon);
                % define the interface in the crust 1.0 model
                m_interface = [];
                for l = 1:size(m0,1)
                    if l == 1
                        m_interface(l,:) = m0(l,:);
                    else
                        m_interface(2*(l-1),:) = [m0(l,1) m0(l-1,2:4)];
                        m_interface(2*(l-1)+1,:) = m0(l,:);
                    end
                end
                dmoho(knode) = m0(end,1);
                
                % Combine with AK135 model below Moho
                dmax = 1000;
                keepz = z > dmoho(knode) & z < dmax;
                m_interface = [m_interface; z(keepz) vp(keepz) vs(keepz) rho(keepz)];
                EPS = 1e-6;
                ztemp = m_interface(:,1);
                idisc = find( ztemp(1:end-1) == ztemp(2:end) );
                ztemp(idisc) = ztemp(idisc) - EPS;
                zpos = 0:dz:zmax;
                vptemp = interp1( ztemp, m_interface(:,2), zpos(1:end-1)+0.5*dz, 'linear','extrap');
                vstemp = interp1( ztemp, m_interface(:,3), zpos(1:end-1)+0.5*dz, 'linear','extrap');
                Z = [Z; zpos(1:end-1)'];
                X = [X; ones(size(vptemp'))*lon];
                Y = [Y; ones(size(vptemp'))*lat];
                VP = [VP; vptemp'];
                VS = [VS; vstemp'];
            end
        end

        % Create interpolants for the 3D velocity model
        Fvp = scatteredInterpolant(X, Y, Z, VP);
        Fvs = scatteredInterpolant(X, Y, Z, VS);

        % Compute velocities at grid points
        [X, Y, Z] = meshgrid(gridStruct.x, gridStruct.y, gridStruct.z);
%         gridStruct.VP = Fvp(,Y,Z);
%         gridStruct.VS = Fvs(X,Y,Z);
        gridStruct.Fvp = Fvp;
        gridStruct.Fvs = Fvs;
end
end