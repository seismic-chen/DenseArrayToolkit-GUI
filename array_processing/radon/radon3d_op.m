function [out] = radon3d_op(in, Param, operator)
% 3D Radon operator (forward/adjoint) for time-domain seismic data
% Inputs:
%   in:     Input data (model or data space)
%   Param:  Structure with parameters:
%           .hx:    x-offset vector [nhx x 1]
%           .hy:    y-offset vector [nhy x 1]
%           .nt:    Number of time samples
%           .dt:    Time sampling interval
%           .type:  Moveout type (1=linear, 2=parabolic, 3=hyperbolic)
%           .px:    x-slopes (linear) / x-curvatures (parabolic) [npx x 1] (types 1-2)
%           .py:    y-slopes (linear) / y-curvatures (parabolic) [npy x 1] (types 1-2)
%           .v:     Velocity vector [nv x 1] (type 3)
%   operator:  1 = forward (model to data), -1 = adjoint (data to model)
% Output:
%   out:    Output data (data or model space)

hx = Param.hx;
hy = Param.hy;
nt = Param.nt;
dt = Param.dt;
type = Param.type;

nhx = length(hx);
nhy = length(hy);

% Initialize based on operator and type
if operator == 1      % Forward operator
    if type == 3      % Hyperbolic
        nv = length(Param.v);
        d = zeros(nt, nhx, nhy);
        m = in;  % Model: [nt, nv]
    else              % Linear/Parabolic
        npx = length(Param.px);
        npy = length(Param.py);
        d = zeros(nt, nhx, nhy);
        m = in;  % Model: [nt, npx, npy]
    end
else                % Adjoint operator
    d = in;  % Data: [nt, nhx, nhy]
    if type == 3
        nv = length(Param.v);
        m = zeros(nt, nv);
    else
        npx = length(Param.px);
        npy = length(Param.py);
        m = zeros(nt, npx, npy);
    end
end

% Perform Radon transform based on type
switch type
    case 1  % Linear (planar surfaces)
        px = Param.px;
        py = Param.py;
        for itau = 1:nt
            for ipx = 1:length(px)
                for ipy = 1:length(py)
                    for ihx = 1:nhx
                        for ihy = 1:nhy
                            t = (itau-1)*dt + px(ipx)*hx(ihx) + py(ipy)*hy(ihy);
                            it = floor(t/dt) + 1;
                            if it >= 1 && it <= nt
                                if operator == 1
                                    d(it, ihx, ihy) = d(it, ihx, ihy) + m(itau, ipx, ipy);
                                else
                                    m(itau, ipx, ipy) = m(itau, ipx, ipy) + d(it, ihx, ihy);
                                end
                            end
                        end
                    end
                end
            end
        end

    case 2  % Parabolic (paraboloid surfaces)
        px = Param.px;
        py = Param.py;
        for itau = 1:nt
            for ipx = 1:length(px)
                for ipy = 1:length(py)
                    for ihx = 1:nhx
                        for ihy = 1:nhy
                            t = (itau-1)*dt + px(ipx)*hx(ihx)^2 + py(ipy)*hy(ihy)^2;
                            it = floor(t/dt) + 1;
                            if it >= 1 && it <= nt
                                if operator == 1
                                    d(it, ihx, ihy) = d(it, ihx, ihy) + m(itau, ipx, ipy);
                                else
                                    m(itau, ipx, ipy) = m(itau, ipx, ipy) + d(it, ihx, ihy);
                                end
                            end
                        end
                    end
                end
            end
        end

    case 3  % Hyperbolic (hyperboloid surfaces)
        v = Param.v;
        for itau = 1:nt
            for iv = 1:length(v)
                for ihx = 1:nhx
                    for ihy = 1:nhy
                        r2 = hx(ihx)^2 + hy(ihy)^2;
                        t = sqrt(((itau-1)*dt)^2 + r2/v(iv)^2);
                        it = floor(t/dt) + 1;
                        if it >= 1 && it <= nt
                            if operator == 1
                                d(it, ihx, ihy) = d(it, ihx, ihy) + m(itau, iv);
                            else
                                m(itau, iv) = m(itau, iv) + d(it, ihx, ihy);
                            end
                        end
                    end
                end
            end
        end
end

% Set output based on operator direction
if operator == 1
    out = d;
else
    out = m;
end
end