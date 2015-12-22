function doit
    R = 0.25;
    R_1 = 0.40;
    delta = 0.2;
    tau = 1e-6;
    h_x = 1/20;
    h_phi = 1/20;
    n_pr = 1.46;
    F_u0 = 1000;
    alpha = 0.01;
    tau_u = 1;

    z_first = R/R_1;
    z_last= 1;

    a = (pi/2 - delta)/(1 - R/R_1);

    y_m = tan(pi/2 - delta);

    function r = z(x)
        r = 1 + 1/a * atan((x - 1) * y_m);
    end

    function r = p(x)
        r = a/y_m .* (1 + (x - 1).^2 * y_m^2);
    end

    function r = lambda(T)
        r = 0.134e-1 * (1 + 4.35e-4*T);
    end

    function r = C(T)
        r = 2.049 + 0.563e-3*T - 0.528e5/T.^2;
    end

    T_ = [293, 1278, 1528, 1677];
    k_p_ = [2e-3, 5e-3, 7.8e-3, 1e-2];
    pp_k = splinefit(T_, k_p_, length(T_), "order", 3);

    function r = F_u(t)
        r = F_u0 * t * exp(-alpha*t / tau_u);
    end

    function r = g(T, t, r)
        sigma = 5.67e-8;
        k_p = ppval(pp_k, T);
        r = F_u(t) * R * k_p / r * exp(-k_p * (R_1 - r)) - 4*k_p * n_pr**2 * sigma * T**4;
    end

    % function r = p_tilde(z, x)
    %     r = ((1 + (x-1)**2 * y_m**2) * atan((x-1)**y_m))/((z-1)*y_m)
    % end

    function A = A_1(T)
        A = zeros(1/h_x, 1/h_phi);
        x = linspace(0, 1, 1/h_x);
        phi = linspace(0, pi/2, 1/h_phi);

        pp_z = splinefit(x, z(x), length(x), "order", 3);
        pp_p = splinefit(x, p(x), length(x), "order", 3);

        L = lambda(T);
        for j = 1:1/h_phi
            for i = 2:1/h_x - 1
                z_half = ppval(pp_z, i - h_x);
                p_half = ppval(pp_p, i - h_x);

                pp_l = splinefit(T(j, :), L(j, :), 1/h_x, "order", 3);
                T_ij = T(j, i);
                T_i1_n = T(j, i - 1);

                l_half_n = ppval(pp_l, linspace(T_i1_n, T_ij, 3)(2));

                A(j, i) = (tau * z_half * p_half * l_half_n) / (h_x * R_1**2);
            end
        end
    end

    function B = B_1(T)
        B = zeros(1/h_x, 1/h_phi);
        x = linspace(0, 1, 1/h_x);
        phi = linspace(0, pi/2, 1/h_phi);

        pp_z = splinefit(x, z(x), length(x), "order", 3);
        pp_p = splinefit(x, p(x), length(x), "order", 3);

        L = lambda(T);

        for j = 1:1/h_phi
            for i = 2:1/h_x - 1
                z_half_n = ppval(pp_z, i - h_x);

                z_half_p = ppval(pp_z, i + h_x);

                p_half_n = ppval(pp_p, i - h_x);
                p_half_p = ppval(pp_p, i + h_x);

                pp_l = splinefit(T(j, :), L(j, :), 1/h_x, "order", 3);

                T_ij = T(j, i);
                T_i1_p = T(j, i + 1);
                T_i1_n = T(j, i - 1);

                l_half_p = ppval(pp_l, linspace(T_ij, T_i1_p, 3)(2));
                l_half_n = ppval(pp_l, linspace(T_i1_n, T_ij, 3)(2));

                c_1 = (z(i) * h_x * C(T(j, i)))/p(i);
                c_2 = (tau * z_half_n * p_half_n * l_half_n)/(h_x * R_1**2);
                c_3 = (tau * z_half_p * p_half_p * l_half_p)/(h_x * R_1**2);

                B(j, i) = c_1 + c_2 + c_3;
            end
        end
    end

    function D = D_1(T)
        C = zeros(1/h_x, 1/h_phi);
        x = linspace(0, 1, 1/h_x);
        phi = linspace(0, pi/2, 1/h_phi);

        pp_z = splinefit(x, z(x), length(x), "order", 3);
        pp_p = splinefit(x, p(x), length(x), "order", 3);

        L = lambda(T);
        for j = 1:1/h_phi
            for i = 2:1/h_x - 1
                z_half = ppval(pp_z, i + h_x);
                p_half = ppval(pp_p, i + h_x);

                pp_l = splinefit(T(j, :), L(j, :), 1/h_x, "order", 3);
                T_ij = T(j, i);
                T_i1_p = T(j, i + 1);

                l_half_p = ppval(pp_l, linspace(T_ij, T_i1_p, 3)(2));

                D(j, i) = (tau * z_half * p_half * l_half_p) / (h_x * R_1**2);
            end
        end
    end

    function F = F_1(T, t, Y)
        F = zeros(1/h_x, 1/h_phi);
        x = linspace(0, 1, 1/h_x);
        phi = linspace(0, pi/2, 1/h_phi);

        % r = linspace(R_1, R, length(x));
        % alpha = linspace(0, pi/2, length(phi));

        for j = 1:1/h_phi
            for i = 2:1/h_x - 1
                z_i = z(i);
                c_ij = C(T(j, i));
                p_i = p(i);

                c_1 = (z_i * h_x * c_ij * Y(j, i))/(p_i);
                c_2 = (tau * z_i * h_x) / p_i * g(T(j, i), t, z_i*R_1);
                F(j, i) = c_1 + c_2;
            end
        end
    end

    function Y = tma(A, B, C, D)
        n = length(B);

        u(1) = C(1) / B(1);
        for i = 2:n-1
            u(i) = C(i)/(B(i) - A(i)*u(i-1));
        end

        v = ones(1, n);

        d(1) = D(1)/B(1);
        for i = 2:n
            d(i) = (D(i) - A(i)*d(i-1))/(B(i) - A(i)*u(i-1));
        end

        Y(n) = d(n);
        for i = fliplr(1:n - 1)
            Y(i) = d(i) - u(i)*Y(i+1);
        end
    end

    function r = tube(t, T)
        A = A_1(T);
        B = B_1(T);
        D = D_1(T);
        Y = F_1(T, t, T);

        for j = 1:1/h_phi
            a_vec = A(j, 2:1/h_x -1);
            b_vec = B(j, 2:1/h_x -1);
            c_vec = D(j, 2:1/h_x -1);
            y_vec = Y(j, 2:1/h_x -1);

            T(j, :) = [293 tma(a_vec, -b_vec, c_vec, -y_vec) 293];
        end
        r = T;
    end

    T = ones(1/h_x, 1/h_phi) .* 293;
    % A_1(T)
    % B_1(T)
    % D_1(T)
    %
    % F_1(T, 0, 0)

    r = tube(tau, T);
    plot(1:columns(r), r(1, :));
end

