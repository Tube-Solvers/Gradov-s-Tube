global R = 0.25;
global R_1 = 0.40;
global delta = 0.2;
global tau = 1e-6;
global h_x = 1/5;
global h_phi = 1/5;
global n_pr = 1.46;

function r = a()
    global delta;
    global R;
    global R_1;

    r = (pi/2 - delta)/(1 - R/R_1);
end

function r = y_m()
    global delta;

    r = tan(pi/2 - delta);
end

function r = z(x)
    r = 1 .+ 1 ./ a .* atan((x .- 1) .* y_m);
end

function r = p(x)
    r = a/y_m .* (1 .+ (x.-1).**2 .* y_m.**2);
end

function r = lambda(T)
    r = 0.134./10 .* (1 .+ 4.35.*1e-4 .* T);
end

function r = C(T)
    r = 2.049 .+ 0.563.*1e-3.*T - 0.528.*1e5./T./T;
end

T_ = [293, 1278, 1528, 1677];
k_p_ = [2e-3, 5e-3, 7.8e-3, 1e-2];
pp = splinefit(T_, k_p_, length(T_), "order", 3);

function r = k_p(T)
    r = ppval(pp, T);
end

function r = g(T)
    sigma = 5.67e-8;
    tau_u = 1;
    F_0u = 1000;
    alpha = 0.01;
    % r = F_0u *
end

function A = A_1(T)
    global tau;
    global h_x;
    global h_phi;
    global R_1;

    A = zeros(1/h_x, 1/h_phi);
    x = linspace(0, 1, 1/h_x);
    phi = linspace(0, pi/2, 1/h_phi);

    pp_z = splinefit(x, z(x), length(x), "order", 3);
    pp_p = splinefit(x, p(x), length(x), "order", 3);

    L = lambda(T);
    for iter_i = 2:1/h_x - 1
        for iter_j = 1:1/h_phi
            z_half = ppval(pp_z, iter_i - h_x);
            p_half = ppval(pp_p, iter_i + h_x);

            pp_l = splinefit(T(iter_j, :), L(iter_j, :), 1/h_x, "order", 3);
            T_ij = T(iter_j, iter_i);
            T_i1_n = T(iter_j, iter_i - 1);

            l_half_n = ppval(pp_l, linspace(T_i1_n, T_ij, 3)(2));

            A(iter_j, iter_i) = (tau * z_half * p_half * l_half_n) / (h_x * R_1**2);
        end
    end
end

function B = B_1(T)
    global tau;
    global h_x;
    global h_phi;
    global R_1;

    B = zeros(1/h_x, 1/h_phi);
    x = linspace(0, 1, 1/h_x);
    phi = linspace(0, pi/2, 1/h_phi);

    pp_z = splinefit(x, z(x), length(x), "order", 3);
    pp_p = splinefit(x, p(x), length(x), "order", 3);

    L = lambda(T);
    for iter_i = 2:1/h_x - 1
        for iter_j = 1:1/h_phi
            z_half_n = ppval(pp_z, iter_i - h_x);
            z_half_p = ppval(pp_z, iter_i + h_x);

            p_half_n = ppval(pp_p, iter_i - h_x);
            p_half_p = ppval(pp_p, iter_i + h_x);

            pp_l = splinefit(T(iter_j, :), L(iter_j, :), 1/h_x, "order", 3);

            T_ij = T(iter_j, iter_i);
            T_i1_p = T(iter_j, iter_i + 1);
            T_i1_n = T(iter_j, iter_i - 1);

            l_half_p = ppval(pp_l, linspace(T_ij, T_i1_p, 3)(2));
            l_half_n = ppval(pp_l, linspace(T_i1_n, T_ij, 3)(2));

            c_1 = (z(iter_i) * h_x * C(T(iter_j, iter_i)))/p(iter_i);
            c_2 = (tau * z_half_n * p_half_n * l_half_p)/(h_x * R_1**2);
            c_3 = (tau * z_half_p * p_half_p * l_half_p)/(h_x * R_1**2);

            B(iter_j, iter_i) = c_1 + c_2 + c_3;
        end
    end
end

function D = D_1(T)
    global tau;
    global h_x;
    global h_phi;
    global R_1;

    C = zeros(1/h_x, 1/h_phi);
    x = linspace(0, 1, 1/h_x);
    phi = linspace(0, pi/2, 1/h_phi);

    pp_z = splinefit(x, z(x), length(x), "order", 3);
    pp_p = splinefit(x, p(x), length(x), "order", 3);

    L = lambda(T);
    for iter_i = 2:1/h_x - 1
        for iter_j = 1:1/h_phi
            z_half = ppval(pp_z, iter_i + h_x);
            p_half = ppval(pp_p, iter_i + h_x);

            pp_l = splinefit(T(iter_j, :), L(iter_j, :), 1/h_x, "order", 3);
            T_ij = T(iter_j, iter_i);
            T_i1_p = T(iter_j, iter_i + 1);

            l_half_p = ppval(pp_l, linspace(T_ij, T_i1_p, 3)(2));

            D(iter_j, iter_i) = (tau * z_half * p_half * l_half_p) / (h_x * R_1**2);
        end
    end
end

T = ones(1/h_x, 1/h_phi) .* 293;
A_1(T)
B_1(T)
D_1(T)
