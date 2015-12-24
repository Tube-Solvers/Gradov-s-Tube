function doit
    R = 0.25;
    R_1 = 0.40;
    delta = 0.2;
    tau = 1e-1;
    n_pr = 1.46;
    F_u0 = 1000;
    alpha = 0.01;
    tau_u = 0.5;
    T_p = 500;
    T_e = 293;

    a = (pi/2 - delta)/(1 - R/R_1);
    y_m = tan(pi/2 - delta);

    k_x = 101;
    h_x = 1/(k_x-1);
    x = (0:k_x-1)*h_x;

    k_phi = 101;
    h_phi = pi/(2*(k_phi-1));
    phi = ((0:k_phi-1)*h_phi)';

    function xl = ml(x)
        xl = [1.5*x(:, 1) - 0.5*x(:, 2), (x(:, 1:end-1) + x(:, 2:end))/2];
    end

    function xr = mr(x)
        xr = [(x(:, 1:end-1) + x(:, 2:end))/2, 1.5*x(:, end) - 0.5*x(:, end-1)];
    end

    xl = ml(x);
    xr = mr(x);

    % Точно ли это z юзается дальше? ещё есть z = r/R_1
    function z = zf(x)
        z = 1 + 1/a * atan((x - 1) * y_m);
    end
    z = zf(x);
    zl = zf(xl);
    zr = zf(xr);

    function p = pf(x)
        p = a/y_m .* (1 + (x - 1).^2 * y_m^2);
    end
    p = pf(x);
    pl = pf(xl);
    pr = pf(xr);

    function r = lambda(T)
        r = 0.134e-1 * (1 + 4.35e-4*T);
    end

    function r = C(T)
        r = 2.049 + 0.563e-3*T - 0.528e5./T.^2;
    end

    function r = F_u(t)
    %    r = F_u0 .* t .* exp(-(alpha.*t).^2./tau_u);
        r = 1000;
    end

    k_t_ = [293, 1278, 1528, 1677];
    k_p_ = [2e-3, 5e-3, 7.8e-3, 1e-2];
    pp_k = splinefit(k_t_, k_p_, length(k_t_), "order", 3);

    function g = g(z, phi, t, T)
        % sigma = 5.67e-8;
        % r = repmat(linspace(R, R_1, columns(z)), rows(phi), 1);
        % k_p = squeeze(ppval(pp_k, T));
        %
        % g = F_u(tau*t) .* R .* k_p ./ r .* exp(-k_p .* (R_1 .- r)) - 4 .* k_p .* n_pr**2 .* sigma .* T.^4;
        % % g = zeros(size(T));
        %
        % g(1:50, :) = 0;

        g = zeros(size(T));
        if t*tau <= 2*tau
            g(25:50, 1) = -1e6;
        end
    end

    function Th = step_x(t, T, Th)
        U = h_x .* z ./ p;

        A = tau./(h_x .* R_1^2) .* zl .* pl .* lambda(ml(Th));
        D = tau./(h_x .* R_1^2) .* zr .* pr .* lambda(mr(Th));
        B = A .+ U.*C(Th) .+ D;
        F = (C(Th).*T .- tau.*g(z, phi, t, Th)) .* U;

        A(:, 1) = nan;
        D(:, 1) = 0;
        B(:, 1) = 1;
        F(:, 1) = Th(:, 1);

        A(:, end) = 0;
        D(:, end) = nan;
        B(:, end) = 1;
        F(:, end) = Th(:, end);

        Th = solve_tridiag(A, -B, D, -F);
    end

    function Th = step_y(t, T, Th)
        V = z .* h_phi;

        A = tau./(z .* R_1^2) .* lambda(ml(Th')');
        D = tau./(z .* R_1^2) .* lambda(mr(Th')');
        B = A .+ V.*C(Th) .+ D;
        F = (C(Th).*T .- tau.*g(z, phi, t, Th)) .* V;

        A(1, :) = nan;
        D(1, :) = 1;
        B(1, :) = 1;
        F(1, :) = 0;

        A(end, :) = 1;
        D(end, :) = nan;
        B(end, :) = 1;
        F(end, :) = 0;

        Th = solve_tridiag(A', -B', D', -F')';
    end

    function Th = iterate_step(t, T, step)
        Th = T;
        while 1
            Th_old = Th;

            Th = step(t, T, Th);

            if max(abs(1 - Th(:)./Th_old(:))) <= 1e-6
                break
            end
        end
    end

    function plot_temperature(i, T)
        printf("Time %d\n", i)
        fig = figure('visible', 'off');
        imagesc(linspace(R, R_1, columns(T)), linspace(pi/2, 0, rows(T)), T, [273, 19000]);
        colorbar();
        print(["plots/" mat2str(t) ".png"], "-dpng")
    end

    function plot_tube(T)
        R_vec = linspace(0.25, 0.40, columns(T));
        P_vec = linspace(0, pi/2, rows(T))';
        X = repmat(R_vec, rows(T), 1) .* cos(repmat(P_vec, 1, columns(T)));
        Y = repmat(R_vec, rows(T), 1) .* sin(repmat(P_vec, 1, columns(T)));

        X = reshape(X, 1, []);
        Y = reshape(Y, 1, []);
        Z = reshape(T, 1, []);

        [xx, yy] = meshgrid(linspace(0, 0.4, 100));
        griddata(X, Y, Z, xx, yy);
    end

    T = T_e * ones(k_phi, k_x);
    % T(:, [1, 2]) = T_p;
    % hold on;
    for t = 1:100
        T = iterate_step(t, T, @(t, T, Th) step_x(t, T, Th));
        T = iterate_step(t, T, @(t, T, Th) step_y(t, T, Th));
        plot_temperature(t, T);
    end

    % plot_tube(T)
end

function x = solve_tridiag(a, b, c, d)
    n = columns(a);

    cc(:, 1) = c(:, 1) ./ b(:, 1);
    for i = 2:n-1
        cc(:, i) = c(:, i)./(b(:, i) - a(:, i).*cc(:, i-1));
    end

    dd(:, 1) = d(:, 1) ./ b(:, 1);
    for i = 2:n
        dd(:, i) = (d(:, i) - a(:, i).*dd(:, i-1)) ./ (b(:, i) - a(:, i).*cc(:, i-1));
    end

    x(:, n) = dd(:, n);
    for i = fliplr(1:n - 1)
        x(:, i) = dd(:, i) - cc(:, i).*x(:, i+1);
    end
end
