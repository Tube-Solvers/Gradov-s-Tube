function doit
    R = 0.25;
    R_1 = 0.40;
    delta = 0.2;
    tau = 1e-2;
    n_pr = 1.46;
    F_u0 = 1000;
    alpha = 0.01;
    tau_u = 0.5;
    T_p = 500;
    T_e = 293;

    a = (pi/2 - delta)/(1 - R/R_1);
    y_m = tan(pi/2 - delta);

    k_x = 21;
    h_x = 1/(k_x-1);
    x = (0:k_x-1)*h_x;

    k_phi = 21;
    h_phi = pi/(2*(k_phi-1));
    phi = ((0:k_phi-1)*h_phi)';

    function xl = ml(x)
        xl = [1.5*x(:, 1) - 0.5*x(:, 2), (x(:, 1:end-1) + x(:, 2:end))/2];
    end

    function xr = mr(x)
        xr = fliplr(ml(fliplr(x)));
    end

    function xu = mu(x)
        xu = ml(x')';
    end

    function xd = md(x)
        xd = mr(x')';
    end

    xl = ml(x);
    xr = mr(x);

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

    %function r = F_u(t)
    %    r = F_u0 .* t .* exp(-(alpha.*t).^2./tau_u);
    %    r = 1000;
    %end

    %k_t_ = [293, 1278, 1528, 1677];
    %k_p_ = [2e-3, 5e-3, 7.8e-3, 1e-2];
    %pp_k = splinefit(k_t_, k_p_, length(k_t_), "order", 3);

    function g = g(t, T)
        % sigma = 5.67e-8;
        % r = repmat(linspace(R, R_1, columns(z)), rows(phi), 1);
        % k_p = squeeze(ppval(pp_k, T));
        %
        % g = F_u(tau*t) .* R .* k_p ./ r .* exp(-k_p .* (R_1 .- r)) - 4 .* k_p .* n_pr**2 .* sigma .* T.^4;
        % % g = zeros(size(T));
        %
        % g(1:50, :) = 0;

        g = zeros(size(T));
        if t < 3*tau
            %g(1:20, 1) = -1e6;
            jm = floor(rows(g) / 2);
            im = floor(columns(g) / 2);
            g(jm-1:jm+1, im-1:im+1) = -1e4;
        end
    end

    function Th = iterate_step(t, T, precision, step)
        Th = T;
        while 1
            Th_old = Th;

            Th = step(t, T, Th);

            if max(abs(1 - Th(:)./Th_old(:))) <= precision
                break
            end
        end
    end

    function Th = step_local_1d_x(t, T, Th)
        U = h_x .* z ./ p;

        A = tau./(h_x .* R_1^2) .* zl .* pl .* lambda(ml(Th));
        D = tau./(h_x .* R_1^2) .* zr .* pr .* lambda(mr(Th));
        B = A .+ U.*C(Th) .+ D;
        F = (C(Th).*T .- tau.*g(t, Th)) .* U;

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

    function Th = step_local_1d_y(t, T, Th)
        V = z .* h_phi;

        A = tau./(z .* R_1^2) .* lambda(mu(Th));
        D = tau./(z .* R_1^2) .* lambda(md(Th));
        B = A .+ V.*C(Th) .+ D;
        F = (C(Th).*T .- tau.*g(t, Th)) .* V;

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

    function Th = step_local_1d(t, T)
        Th = T;

        Th = iterate_step(t, Th, 1e-6, @(t, T, Th) step_local_1d_x(t, T, Th));
        Th = iterate_step(t, Th, 1e-6, @(t, T, Th) step_local_1d_y(t, T, Th));
    end

    function Th = step_naive_2d_xy(t, T, Th)
        H = x.*h_x.*h_phi./p;

        A = tau.*h_x./(R_1^2 .* p .* z) .* lambda(mu(Th));
        B = tau.*h_phi./(R_1^2 .* h_x) .* zl .* pl .* lambda(ml(Th));
        D = tau.*h_phi./(R_1^2 .* h_x) .* zr .* pr .* lambda(mr(Th));
        E = tau.*h_x./(R_1^2 .* p .* z) .* lambda(md(Th));
        G = A .+ B .+ H.*C(Th) .+ D .+ E;
        F = (C(Th).*T - tau.*g(t, Th)) .* H;

        A(1, :) = nan;
        A(end, :) = 1;

        B(:, 1) = nan;
        B(:, end) = 0;

        D(:, 1) = 0;
        D(:, end) = nan;

        E(1, :) = 1;
        E(end, :) = nan;

        G([1, end], :) = 1;
        G(:, [1, end]) = 1;

        F([1, end], :) = 0;
        F(:, [1, end]) = Th(:, [1, end]);
        
        L = [];
        R = [];
        for j = 1:k_phi
            for i = 1:k_x
                if j > 1
                    L(i, j, i, j-1) = A(j, i);
                end
                if i > 1
                    L(i, j, i-1, j) = B(j, i);
                end
                L(i, j, i, j) = -G(j, i);
                if i < k_x
                    L(i, j, i+1, j) = D(j, i);
                end
                if j < k_phi
                    L(i, j, i, j+1) = E(j, i);
                end
                R(i, j) = -F(j, i);
            end
        end

        L = reshape(L, k_phi * k_x, k_phi * k_x);
        R = reshape(R, k_phi * k_x, 1);

        Th = L \ R;

        Th = reshape(Th, k_phi, k_x);
    end

    function Th = step_naive_2d(t, T)
        Th = iterate_step(t, T, 1e-2, @(t, T, Th) step_naive_2d_xy(t, T, Th));
    end

    function plot_temperature(t, T, file, delay)
        graphics_toolkit('gnuplot');
        printf('Time %d\n', t/tau);
        fflush(stdout);
        figure('visible', 'off');
        lo = min(T(:));
        hi = max(T(:)); hi = max(hi, lo + 1);
        imagesc(linspace(R, R_1, columns(T)), linspace(0, 2*pi*R, rows(T)), T, [lo, hi]);
        colorbar();
        print('.temp.png', '-dpng');
        [data, map] = rgb2ind(imread('.temp.png'));
        unlink('.temp.png');
        if t == 0
            imwrite(data, map, file, 'gif', 'LoopCount', Inf,
                                            'DelayTime', delay); 
        else
            imwrite(data, map, file, 'gif', 'WriteMode', 'append',
                                            'DelayTime', delay);
        end
    end

    T = T_e * ones(k_phi, k_x);
    for t = (0:99).*tau
        T = step_local_1d(t, T);
        %T = step_naive_2d(t, T);
        plot_temperature(t, T, 'animation.gif', 0.1);
    end
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
