clear;clc;clf;close all;

run('quarter_80.m')

miu  =  0.3; % Poisson's ratio
E    = 1E9; % Elastic Modulus

% D matrix
D = zeros(3, 3);                         
D(1, 1) = E/(1-miu^2);
D(2, 2) = D(1,1);
D(1, 2) = E *miu/(1-miu^2);
D(2, 1) = D(1,2);
D(3, 3) = E*(1-miu)/(2*(1-miu^2));

% 2D quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% 1D quadrature rule
n_g = 5;
[seit, seit_weight] = Gauss(n_g , -1 ,1);

% mesh generation
n_en   = 4;                    % number of nodes in an element
n_el   = length(msh.QUADS);    % total number of elements
n_np   = msh.nbNod;            % total number of nodal points
n_lin = length(msh.LINES);     % total number of lines

% generate the nodal coordinates
x_coor = msh.POS(:,1);         % x-axis POS of all nodal points
y_coor = msh.POS(:,2);         % y-axis POS of all nodal points

% IEN array
IEN_tri = zeros(1,1);
IEN = msh.QUADS(:,1:4);

for ee = 1:size(IEN,1)
    IEN_tri(ee*2-1,1) = IEN(ee,1);
    IEN_tri(ee*2-1,2) = IEN(ee,2);
    IEN_tri(ee*2-1,3) = IEN(ee,3);

    IEN_tri(ee*2,1) = IEN(ee,1);
    IEN_tri(ee*2,2) = IEN(ee,3);
    IEN_tri(ee*2,3) = IEN(ee,4);
end

% 右边和上边条件设置
Dirichlet = -1;
Newman    = -1;

% ID array
ID = -1 * ones(n_np, 2);
for ii = 1 : n_lin
    if msh.LINES(ii , 3) == 8                 % 右边界
        ID(msh.LINES(ii,1),1) = Dirichlet;
        ID(msh.LINES(ii,1),2) = Dirichlet;
    elseif msh.LINES(ii , 3) == 9             % 上边界
        ID(msh.LINES(ii,1),1) = Newman;
        ID(msh.LINES(ii,1),2) = Newman;
    end

    if msh.LINES(ii , 3) == 10                % 左边界
        ID(msh.LINES(ii,1),1) = 0;
    elseif msh.LINES(ii , 3) == 11            % 下边界
        ID(msh.LINES(ii,1),2) = 0;
    end
end

ID(3,1) = 0;
ID(4,2) = 0;

counter = 1;
for ii = 1 : n_np
    for jj = 1 : 2
        if ID(ii,jj) ~= 0
            ID(ii,jj) = counter;
            counter = counter + 1;
        end
    end
end

n_eq = counter -1;

LM = cell(2, 1);                         % need to use a new map for LM
for ii = 1 : 2
    LM{ii} = zeros(size(IEN));           % so use cell
    for jj = 1 : n_el
        for kk = 1 : n_en
            cc = IEN(jj, kk);
            LM{ii}(jj, kk)= ID(cc, ii);
        end
    end
end

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
% in order to conveniently use coordinates to determine whether the unit is
% on the boundary, let x_ele(1, n_en) change into (x1,x2,x3,x4,x1)
x_ele = zeros(1, n_en + 1);
y_ele = zeros(1, n_en + 1);

for ee = 1 : n_el
    x_ele(1:n_en) = x_coor( IEN(ee, 1:n_en) ); x_ele(5) = x_ele(1);
    y_ele(1:n_en) = y_coor( IEN(ee, 1:n_en) ); y_ele(5) = y_ele(1);

    k_ele = zeros(n_en * 2, n_en * 2); % element stiffness matrix
    f_ele = zeros(n_en * 2, 1);    % element load vector

    for ii = 1 : n_g
        y_l = 0.0; dy_dseit = 0.0;
        x_l = 0.0; dx_dseit = 0.0;
        for aa = 1 : n_en
            pp = 2 * aa;
            if x_ele(aa) == 1 && x_ele(aa + 1) == 1
                for jj = 1 : 2
                    y_l      = y_l      + y_ele(aa + jj - 1) * PolyShape(1, jj , seit(ii),0);
                    dy_dseit = dy_dseit + y_ele(aa + jj - 1) * PolyShape(1, jj , seit(ii),1);
                end
                hx = 10000; hy = 0;
                for jj = 1 : 2
                    f_ele(pp - 1) = f_ele(pp - 1) + seit_weight(ii) * PolyShape(1, jj , seit(ii),0) * hx * dy_dseit;
                    f_ele(pp) = f_ele(pp) + seit_weight(ii) * PolyShape(1, jj , seit(ii),0) * hy * dy_dseit;
                end
            end
            % 
            % if y_ele(aa) == 1 && y_ele(aa + 1) == 1
            %     for jj = 1 : 2
            %         x_l      = x_l      + x_ele(aa + jj - 1) * PolyShape(1, jj , seit(ii),0);
            %         dx_dseit = dx_dseit + x_ele(aa + jj - 1) * PolyShape(1, jj , seit(ii),1);
            %     end
            %     [hx , hy] = fromPointToHi(x_l,1,0,1);
            %     for jj = 1 : 2
            %         f_ele(pp - 1) = f_ele(pp - 1) + seit_weight(ii) * PolyShape(1, jj , seit(ii),0) * hx * dy_dseit;
            %         f_ele(pp) = f_ele(pp) + seit_weight(ii) * PolyShape(1, jj , seit(ii),0) * hy * dy_dseit;
            %     end
            % end
        end
    end

    for aa = 1 : n_en
        for ii = 1 : 2
            PP = LM{ii}(ee,aa);
            if PP > 0
                F(PP) = F(PP) + f_ele(2*(aa-1)+ii);
            end
        end
    end

    for ll = 1 : n_int
        x_l = 0.0; y_l = 0.0;
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        for aa = 1 : n_en
            x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
            y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end

        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;


        for aa = 1 : n_en
            Na = Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

            B1 = zeros(3, 2);B1(1, 1) = Na_x;B1(2, 2) = Na_y;B1(3, 2) = Na_x;B1(3, 1) = Na_y;

            for bb = 1 : n_en
                Nb = Quad(bb, xi(ll), eta(ll));
                [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;

                B2 = zeros(3, 2);B2(1, 1) = Nb_x;B2(2, 2) = Nb_y;B2(3, 2) = Nb_x;B2(3, 1) = Nb_y;

                for i = 1 : 2
                    for j = 1 : 2
                        ei = unit_vector(i);
                        ej = unit_vector(j);
                        k_ele(2*(aa-1)+i, 2*(bb-1)+j) = k_ele(2*(aa-1)+i, 2*(bb-1)+j) + weight(ll) * detJ * ei' * B1' * D * B2 * ej;
                    end
                end
            end % end of bb loop
        end % end of aa loop
    end % end of quadrature loop

    for aa = 1 : n_en
        for ii = 1 : 2
            PP = LM{ii}(ee,aa);
            if PP > 0
                % F(PP) = F(PP) + f_ele(2*(aa-1)+ii);
                for bb = 1 :n_en
                    for jj = 1: 2
                        QQ = LM{jj}(ee,bb);
                        if QQ > 0
                            K(PP, QQ) = K(PP, QQ) + k_ele(2*(aa-1)+ii, 2*(bb-1)+jj);
                        else

                        end
                    end
                end
            end
        end
    end
end

dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, 2);
for n_nsd = 1 : 2
    for ii = 1 : n_np
        index = ID(ii , n_nsd);
        if index > 0
            disp(ii, n_nsd) = dn(index);
        else

        end
    end
end

trisurf(IEN_tri, x_coor+1, y_coor+1, disp(:,1));
axis equal;
colormap jet
shading interp
title("x方向位移")
view(2);

figure;
trisurf(IEN_tri, x_coor+1, y_coor+1, disp(:,2));
axis equal;
colormap jet
shading interp
title("y方向位移")
view(2);


% another quadrature rule
xi  = [ -1 1 1 -1];
eta = [ -1 -1 1 1];
dux_dx = zeros(n_np , 1); duy_dy = dux_dx; dux_dy = dux_dx; duy_dx = dux_dx;

for ll = 1 :n_en                   % 类高斯积分点
    for ee = 1 : n_el              % 单元点
        x_ele = x_coor( IEN(ee, 1:n_en) );
        y_ele = y_coor( IEN(ee, 1:n_en) );

        u_ele = disp( IEN(ee, :) , 1);
        v_ele = disp( IEN(ee, :) , 2);

        x_l = 0.0; dx_dxi = 0.0; dx_deta = 0.0;
        y_l = 0.0; dy_dxi = 0.0; dy_deta = 0.0;
        du_dx = 0.0; du_dy = 0; dv_dx = 0;dv_dy = 0;
        for aa = 1 : n_en
            x_l     = x_l     + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
            y_l     = y_l     + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end

        detJ   = dx_dxi * dy_deta - dx_deta * dy_dxi;

        for aa = 1 : n_en
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            Na = Quad(aa, xi(ll), eta(ll));

            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
            du_dx = du_dx + Na_x * u_ele(aa);
            du_dy = du_dy + Na_y * u_ele(aa);
            dv_dx = dv_dx + Na_x * v_ele(aa);
            dv_dy = dv_dy + Na_y * v_ele(aa);

        end
        dux_dx(ee,ll) = du_dx; dux_dy(ee,ll) = du_dy;
        duy_dx(ee,ll) = dv_dx; duy_dy(ee,ll) = dv_dy;
    end
end

strain = zeros(n_np, 3); number = zeros(n_np, 1);
for ee = 1 : n_el
    for ll = 1 : n_en
        index = IEN(ee,ll);
        strain(index,1) = strain(index,1) + dux_dx(ee,ll);
        strain(index,2) = strain(index,2) + duy_dy(ee,ll);
        strain(index,3) = strain(index,3) + dux_dy(ee,ll) + duy_dx(ee,ll);
        number(index) = number(index) + 1;
    end
end
for ii = 1 : 3
    strain(:,ii) = strain(:,ii) ./ number;
end
stress = D*strain'; 

for ii = 1 : 3
    figure;
    trisurf(IEN_tri, x_coor+1, y_coor+1, strain(:,ii));
    axis equal;
    colormap jet
    shading interp
    if ii == 1
        title("x方向应变")
    elseif ii == 2
        title("y方向应变")
    elseif ii == 3
        title("xy方向应变")
    end
    view(2);

    figure;
    trisurf(IEN_tri, x_coor+1, y_coor+1, stress(ii,:));
    axis equal;
    colormap jet
    shading interp
    if ii == 1
        title("x方向应力")
    elseif ii == 2
        title("y方向应力")
    elseif ii == 3
        title("xy方向应力")
    end
    view(2);
end

function ei = unit_vector(ii)
if ii == 1
    ei = [1, 0]';
else
    ei = [0, 1]';
end
end