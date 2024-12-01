clear all; clc; clf; % clean the memory, screen, and figure

% Problem definition
f = @(x) -20*x.^3; % f(x) is the source
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0

% Setup the mesh
pp   = 1;              % polynomial degree
n_en = pp + 1;         % number of element or local nodes
n_int = 10;
[xi, weight] = Gauss(n_int, -1, 1);
[t01, t01_weight] = Gauss(n_int, 0, 1);
eL21 = zeros(8, 1);
eH11 = zeros(8, 1);

for n_el = 2 : 2 : 16                             % hw4 2
    n_np = n_el * pp + 1;                         % copy above which changed by n_el
    n_eq = n_np - 1;

    hh = 1.0 / (n_np - 1);
    x_coor = 0 : hh : 1;

    IEN = zeros(n_el , n_en);
    for ee = 1 : n_el
        for aa = 1 : n_en
            IEN(ee , aa) = (ee - 1) * pp + aa;
        end
    end

    ID = 1 : n_np;
    ID(end) = 0;


    K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq);
    F = zeros(n_eq, 1);

    for ee = 1 : n_el
        k_ele = zeros(n_en, n_en);
        f_ele = zeros(n_en, 1);

        x_ele = x_coor(IEN(ee,:));

        for qua = 1 : n_int    
            dx_dxi = 0.0;
            x_l = 0.0;
            for aa = 1 : n_en
                x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0);
                dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
            end
            dxi_dx = 1.0 / dx_dxi;

            for aa = 1 : n_en
                f_ele(aa) = f_ele(aa) + weight(qua) * PolyShape(pp, aa, xi(qua), 0) * f(x_l) * dx_dxi;
                for bb = 1 : n_en
                    k_ele(aa, bb) = k_ele(aa, bb) + weight(qua) * PolyShape(pp, aa, xi(qua), 1) * PolyShape(pp, bb, xi(qua), 1) * dxi_dx;
                end
            end
        end

        for aa = 1 : n_en
            P = ID(IEN(ee,aa));
            if(P > 0)
                F(P) = F(P) + f_ele(aa);
                for bb = 1 : n_en
                    Q = ID(IEN(ee,bb));
                    if(Q > 0)
                        K(P, Q) = K(P, Q) + k_ele(aa, bb);
                    else
                        F(P) = F(P) - k_ele(aa, bb) * g; % handles the Dirichlet boundary data
                    end
                end
            end
        end
    end

    F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h;


    d_temp = K \ F;

    disp = [d_temp; g];

    n_sam = 20;
    xi_sam = -1 : (2/n_sam) : 1;

    x_sam = zeros(n_el * n_sam + 1, 1);
    y_sam = x_sam; % store the exact solution value at sampling points
    u_sam = x_sam; % store the numerical solution value at sampling pts

    for ee = 1 : n_el
        x_ele = x_coor( IEN(ee, :) );
        u_ele = disp( IEN(ee, :) );

        if ee == n_el
            n_sam_end = n_sam+1;
        else
            n_sam_end = n_sam;
        end

        for ll = 1 : n_sam_end
            x_l = 0.0;
            u_l = 0.0;
            for aa = 1 : n_en
                x_l = x_l + x_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
                u_l = u_l + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
            end

            x_sam( (ee-1)*n_sam + ll ) = x_l;
            u_sam( (ee-1)*n_sam + ll ) = u_l;
            y_sam( (ee-1)*n_sam + ll ) = x_l^5;
        end
    end


    % from here begin
    for nel = 1 : n_el
        [xi_nel, weighti_nel] = Gauss(n_int, (nel - 1) /n_el , nel / n_el);     %每个节点内部的积分要在哪里，节点左右由间距表示
        x_ele = x_coor( IEN(nel, :) );
        u_ele = disp( IEN(nel, :) );
        for qua = 1 : n_int
            integ = 0.0;
            for aa = 1 : n_en
                integ = integ + u_ele(aa) * PolyShape(pp, aa, xi_nel(qua), 0);
            end
            eL21(n_el/2, 1) = eL21(n_el/2, 1) + weighti_nel(qua) * (integ - xi_nel(qua)^5)^2;
        end
    end

    t01_i = zeros(n_int , 1);
    for ee = 1 : n_el
        for qua = 1 : n_int
            if (1 / n_el) * (ee - 1) < t01(qua) && t01(qua) < (1 / n_el) * ee
                t01_i(qua) = (t01(qua) - (1 / n_el) * (ee - 1)) / (1 / n_el) * 2 - 1;
            end
        end
    end
end


du_dx = @(x) 5 * x.^4;
uu    = @(x) x.^5;
eL22 = 0.0;
eH12 = 0.0;

for qua = 1 : n_int
    eH12 = eH12 + t01_weight(qua) * du_dx(t01(qua)) ^ 2;
    eL22 = eL22 + t01_weight(qua) * uu(t01(qua)) ^ 2;
    if qua == n_int
        eH12 = eH12 ^ 0.5;
        eL22 = eL22 ^ 0.5;
    end
end

eL2 = eL21 / eL22;
eH1 = eH11 / eH12;





























