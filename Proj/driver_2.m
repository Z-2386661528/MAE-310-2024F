clear all; clc;

run('quarter.m')

miu  =  0.3; % Poisson's ratio
E    = 1E9; % Elastic Modulus

D = zeros(3, 3);                         % another D
D(1, 1) = E/(1-miu^2);
D(2, 2) = D(1,1);
D(1, 2) = E *miu/(1-miu^2);
D(2, 1) = D(1,2);
D(3, 3) = E*(1-miu)/(2*(1-miu^2));

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

n_g = 5;
[seit, seit_weight] = Gauss(n_g , -1 ,1);


n_lin = size(msh.LINES);
n_lin = n_lin(1,1);
x_coor = zeros(msh.nbNod, 1);
y_coor = x_coor;

for i = 1 : msh.nbNod
    x_coor(i) = msh.POS(i,1);
    y_coor(i) = msh.POS(i,2);
end


% mesh generation
n_en   = 4;               % number of nodes in an element
n_el   = 98;              % total number of elements
n_np   = msh.nbNod;       % total number of nodal points
n_nd   = 7;

% generate the nodal coordinates
xcoor = msh.POS(:,1);
ycoor = msh.POS(:,2);

% IEN array
%IEN array
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

for i = 1 : n_el/2
    a = IEN(i, 1);
    b = IEN(i, 2);
    IEN(i, 1) = IEN(i, 4);
    IEN(i, 2) = IEN(i, 3);
    IEN(i, 4) = a;
    IEN(i, 3) = b;
end

% 右边和上边条件设置
Dirichlet = -1;
Newman    = -1;

% ID array
ID = -1 * ones(n_np, 2);
for ii = 1 : n_lin
    if line(ii , 3) == 8
        ID(line(ii,1),1) = Dirichlet;
        ID(line(ii,1),2) = Dirichlet;
    elseif line(ii , 3) == 9
        ID(line(ii,1),1) = Newman;
        ID(line(ii,1),2) = Newman;
    end
    if line(ii , 3) == 10
        ID(line(ii,1),1) = 0;
    elseif line(ii , 3) == 11
        ID(line(ii,1),2) = 0;
    end
end

ID(1,1) = Dirichlet; ID(1,2) = 0;
ID(2,1) = Dirichlet; ID(2,2) = Newman;
ID(3,1) = 0;         ID(3,2) = Newman;
ID(4,1) = Newman;    ID(4,2) = 0;
ID(5,1) = 0;         ID(5,2) = Dirichlet;

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

% find the up and right lines
x1counter = 0; nodel_up = zeros(n_nd,3);
y1counter = 0; nodel_right = zeros(n_nd,3);
x2counter = 0; nodel_down = zeros(n_nd,3);
y2counter = 0; nodel_left = zeros(n_nd,3);

for ii = 1 : n_lin
    if msh.LINES(ii,3) == 9
        x1counter = x1counter + 1;
        nodel_up(x1counter,:) = msh.LINES(ii,:);
    elseif msh.LINES(ii,3) == 8
        y1counter = y1counter + 1;
        nodel_right(y1counter,:) = msh.LINES(ii,:);
    elseif msh.LINES(ii,3) == 11
        x2counter = x2counter + 1;
        nodel_down(x2counter,:) = msh.LINES(ii,:);
    elseif msh.LINES(ii,3) == 10
        y2counter = y2counter + 1;
        nodel_left(y2counter,:) = msh.LINES(ii,:);
    end
end

% allocate the stiffness matrix and load vector
% K = spalloc(n_eq, n_eq, 9 * n_eq);
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
x_ele = zeros(n_el, 5);
y_ele = zeros(n_el, 5);
for ee = 1 : n_el
    x_ele(ee,1:n_en) = x_coor( IEN(ee, 1:n_en) );
    x_ele(ee,5) = x_ele(ee,1);
    y_ele(ee,1:n_en) = y_coor( IEN(ee, 1:n_en) );
    y_ele(ee,5) = y_ele(ee,1);

    k_ele = zeros(n_en * 2, n_en * 2); % element stiffness matrix
    f_ele = zeros(n_en * 2, 1);    % element load vector

    for ii = 1 : n_g
        y_l = 0.0; dy_dseit = 0.0;
        x_l = 0.0; dx_dseit = 0.0;
        for aa = 1 : n_en
            pp = 2 * aa - 2;
            if x_ele(aa) == 1 && x_ele(aa + 1) == 1
                for jj = 1 : 2
                    y_l      = y_l      + y_ele(aa + jj - 1) * PolyShape(1, jj , seit(ii),0);
                    dy_dseit = dy_dseit + y_ele(aa + jj - 1) * PolyShape(1, jj , seit(ii),1);
                end

                [hx , hy] = fromPointToHi(1,y_l, 1, 0);
                for jj = 1 : 2
                    f_ele(pp + 1) = f_ele(pp + 1) + seit_weight(ii) * PolyShape(1, jj , seit(ii),0) * hx * dy_dseit;
                    f_ele(pp + 2) = f_ele(pp + 2) + seit_weight(ii) * PolyShape(1, jj , seit(ii),0) * hy * dy_dseit;
                end
            end

            if x_ele(aa) == -1 && x_ele(aa + 1) == -1
                for jj = 1 : 2
                    y_l      = y_l      + y_ele(aa + jj - 1) * PolyShape(1, jj , seit(ii),0);
                    dy_dseit = dy_dseit + y_ele(aa + jj - 1) * PolyShape(1, jj , seit(ii),1);
                end

                [hx , hy] = fromPointToHi(1 , y_l, -1, 0);
                for jj = 1 : 2
                    f_ele(pp + 1) = f_ele(pp + 1) + seit_weight(ii) * PolyShape(1, jj , seit(ii),0) * hx * dy_dseit;
                    f_ele(pp + 2) = f_ele(pp + 2) + seit_weight(ii) * PolyShape(1, jj , seit(ii),0) * hy * dy_dseit;
                end
            end

            if y_ele(aa) == 1 && y_ele(aa + 1) == 1
                for jj = 1 : 2
                    x_l      = x_l      + x_ele(aa + jj - 1) * PolyShape(1, jj , seit(ii),0);
                    dx_dseit = dx_dseit + x_ele(aa + jj - 1) * PolyShape(1, jj , seit(ii),1);
                end

                [hx , hy] = fromPointToHi(1 , y_l, 0, 1);
                for jj = 1 : 2
                    f_ele(pp + 1) = f_ele(pp + 1) + seit_weight(ii) * PolyShape(1, jj , seit(ii),0) * hx * dy_dseit;
                    f_ele(pp + 2) = f_ele(pp + 2) + seit_weight(ii) * PolyShape(1, jj , seit(ii),0) * hy * dy_dseit;
                end
            end

            if y_ele(aa) == -1 && y_ele(aa + 1) == -1
                for jj = 1 : 2
                    x_l      = x_l      + x_ele(aa + jj - 1) * PolyShape(1, jj , seit(ii),0);
                    dx_dseit = dx_dseit + x_ele(aa + jj - 1) * PolyShape(1, jj , seit(ii),1);
                end

                [hx , hy] = fromPointToHi(1 , y_l, 0, -1);
                for jj = 1 : 2
                    f_ele(pp + 1) = f_ele(pp + 1) + seit_weight(ii) * PolyShape(1, jj , seit(ii),0) * hx * dy_dseit;
                    f_ele(pp + 2) = f_ele(pp + 2) + seit_weight(ii) * PolyShape(1, jj , seit(ii),0) * hy * dy_dseit;
                end
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
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dseit) / detJ;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dseit) / detJ;

            B1 = zeros(3, 2);B1(1, 1) = Na_x;B1(2, 2) = Na_y;B1(3, 2) = Na_x;B1(3, 1) = Na_y;


            for bb = 1 : n_en
                Nb = Quad(bb, xi(ll), eta(ll));
                [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dseit) / detJ;
                Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dseit) / detJ;
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
                F(PP) = F(PP) + f_ele(2*(aa-1)+ii);
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

hold on;
trisurf(IEN_tri, x_coor, y_coor, disp(:,1));
axis equal;
colormap jet
shading interp





function ei = unit_vector(ii)
if ii == 1
    ei = [1, 0]';
else
    ei = [0, 1]';
end
end