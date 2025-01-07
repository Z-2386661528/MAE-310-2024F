clear all; clc;

kappa = 1.0; % conductivity

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) 2.0 * kappa * x * (1-x) + 2.0 * kappa  *y * (1-y); % source term

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta, 3);


for i = 1 : 6
    % mesh generation
    n_en   = 3;                    % number of nodes in an element
    n_el_x = 10 * i;               % number of elements in y-dir
    n_el_y = n_el_x;               % number of elements in y-dir
    n_el   = 2 * n_el_x * n_el_y;  % total number of elements

    n_np_x = n_el_x + 1;      % number of nodal points in x-dir
    n_np_y = n_el_y + 1;      % number of nodal points in y-dir
    n_np   = n_np_x * n_np_y; % total number of nodal points

    x_coor = zeros(n_np, 1);
    y_coor = x_coor;

    hx = 1.0 / n_el_x;        % mesh size in x-dir
    hy = 1.0 / n_el_y;        % mesh size in y-dir

    % generate the nodal coordinates
    for ny = 1 : n_np_y
    for nx = 1 : n_np_x
        index = (ny-1)*n_np_x + nx; % nodal index
        x_coor(index) = (nx-1) * hx;
        y_coor(index) = (ny-1) * hy;
    end
    end
 
    % IEN array
    IEN = zeros(n_el, n_en);
    for ex = 1 : n_el_x
        for ey = 1 : n_el_y
            ee = (ey-1) * n_el_x + ex; % element index
    
            IEN(2 * ee - 1, 1) = (ey-1) * n_np_x + ex;         %下直角点1
            IEN(2 * ee - 1, 2) = (ey-1) * n_np_x + ex + 1;     %下直角点2
            IEN(2 * ee - 1, 3) =  ey    * n_np_x + ex;         %下直角点3

            IEN(2 * ee, 1) =   ey * n_np_x + ex + 1;           %上直角点1
            IEN(2 * ee, 2) = IEN(2 * ee - 1, 3);               %下直角点2
            IEN(2 * ee, 3) = IEN(2 * ee - 1, 2);               %下直角点3
        end
    end

    % ID array
    ID = zeros(n_np,1);
    counter = 0;
    for ny = 2 : n_np_y - 1
        for nx = 2 : n_np_x - 1                              %从左往右，再从下到上，所以是先y后x，只需要2-5
        index = (ny-1)*n_np_x + nx;
        counter = counter + 1;
        ID(index) = counter;  
        end
    end

    n_eq = counter;

    LM = ID(IEN);

    % allocate the stiffness matrix and load vector
    K = spalloc(n_eq, n_eq, 9 * n_eq);
    F = zeros(n_eq, 1);

    % loop over element to assembly the matrix and vector
    for ee = 1 : n_el
        x_ele = x_coor( IEN(ee, 1:n_en) );
        y_ele = y_coor( IEN(ee, 1:n_en) );
  
        k_ele = zeros(n_en, n_en); % element stiffness matrix
        f_ele = zeros(n_en, 1);    % element load vector
  
        for ll = 1 : 3         %ll是积分过程中的节点
            x_l = 0.0; y_l = 0.0;
            dx_dxi = 0.0; dx_deta = 0.0;
            dy_dxi = 0.0; dy_deta = 0.0;
                for aa = 1 : n_en
                    x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll), n_en);
                    y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll), n_en);    
                    [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll), n_en);
                    dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
                    dx_deta = dx_deta + x_ele(aa) * Na_eta;
                    dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
                    dy_deta = dy_deta + y_ele(aa) * Na_eta;
                end
    
                detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;              %check
    
            for aa = 1 : n_en
                Na = Quad(aa, xi(ll), eta(ll), n_en);
                [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll), n_en);
                Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
      
                f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;
      
                for bb = 1 : n_en
                    Nb = Quad(bb, xi(ll), eta(ll), n_en);
                    [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll), n_en);
                    Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                    Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
        
                    k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
                end % end of bb loop
            end % end of aa loop
        end % end of quadrature loop
 
        for aa = 1 : n_en
            PP = LM(ee, aa);
            if PP > 0
                F(PP) = F(PP) + f_ele(aa);
      
                for bb = 1 : n_en
                    QQ = LM(ee, bb);
                    if QQ > 0
                        K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
                    else
                        % modify F with the boundary data
                        % here we do nothing because the boundary data g is zero or
                        % homogeneous
                    end
                end  
            end
        end
    end

    % solve the stiffness matrix
    dn = K \ F;

    % insert dn back into the vector for all nodes
    disp = zeros(n_np, 1);

    for ii = 1 : n_np
        index = ID(ii);
        if index > 0
            disp(ii) = dn(index);
        else
            % modify disp with the g data. Here it does nothing because g is zero
        end
    end

    % save the solution vector and number of elements to disp with name
    % HEAT.mat
    save("HEAT", "disp", "n_el_x", "n_el_y");
    [xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta, 3);
    e0 = 0.0; e1 = 0.0;

    for ee = 1 : n_el
        x_ele = x_coor( IEN(ee, :) );
        y_ele = y_coor( IEN(ee, :) );
        u_ele = disp( IEN(ee, :) );

        for ll = 1 : 3
            x_l = 0.0; dx_dxi = 0.0; dx_deta = 0.0;
            y_l = 0.0; dy_dxi = 0.0; dy_deta = 0.0;
            u_l = 0.0; du_dx = 0.0; du_dy = 0.0;
                for aa = 1 : n_en
                    x_l     = x_l     + x_ele(aa) * Quad(aa, xi(ll), eta(ll), n_en);
                    y_l     = y_l     + y_ele(aa) * Quad(aa, xi(ll), eta(ll), n_en);
                    [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll), n_en);
                    dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
                    dx_deta = dx_deta + x_ele(aa) * Na_eta;
                    dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
                    dy_deta = dy_deta + y_ele(aa) * Na_eta;
                end

                detJ   = dx_dxi * dy_deta - dx_deta * dy_dxi;

                for aa = 1 : n_en
                    [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll), n_en);
                    Na = Quad(aa, xi(ll), eta(ll), n_en);
                    Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                    Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
                    u_l = u_l + u_ele(aa) * Na;
                    du_dx = du_dx + u_ele(aa) * Na_x;
                    du_dy = du_dy + u_ele(aa) * Na_y;
                end

                e0 = e0 + weight(ll) * (u_l - exact(x_l, y_l) )^2 * detJ;
                e1 = e1 + weight(ll) * ((du_dx - exact_x(x_l, y_l))^2 + (du_dy - exact_y(x_l, y_l))^2) * detJ;
        end
    end
    e0r(i) = e0 ^ 0.5;
    e1r(i) = e1 ^ 0.5;
end

plot(log(1./(10 : 10 : 60)), log(e0r), '-r','LineWidth',3);
hold on;
grid on;
e00Slope = polyfit(log(1./(10 : 10 : 60)), log(e0r), 1);
e10Slope = polyfit(log(1./(10 : 10 : 60)), log(e1r), 1);

plot(log(1./(10 : 10 : 60)), log(e1r), '-b','LineWidth',3);

str1 = 'e0r  k='+ "" + mat2str(e00Slope(1, 1));
str2 = 'e1r  k='+ "" + mat2str(e10Slope(1, 1));

legend(str1,str2);
title("triangle");

% EOF