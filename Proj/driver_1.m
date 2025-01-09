clear all; clc;

miu  =  0.3; % Poisson's ratio
E    = 1E9; % Elastic Modulus

D = zeros(3, 3);                         % another D  
D(1, 1) = E/(1-miu^2);
D(2, 2) = D(1,1);
D(1, 2) = E *miu/(1-miu^2);
D(2, 1) = D(1,2);
D(3, 3) = E*(1-miu)/(2*(1-miu^2));

% exact solution
exact_at_x = @(x,y) x*(1-x)*y*(1-y);
exact_at_y = @(x,y) x*(1-x)*y*(1-y);

exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f_1 = @(x,y) (2*E*y*(y - 1))/(miu^2 - 1) - (E*(miu/2 - 1/2)*((x - 1)*(y - 1) ...
    + x*y + 2*x*(x - 1) + x*(y - 1) + y*(x - 1)))/(miu^2 - 1) ...
+ (E*miu*((x - 1)*(y - 1) + x*y + x*(y - 1) + y*(x - 1)))/(miu^2 - 1);
f_2 = @(x,y) (2*E*x*(x - 1))/(miu^2 - 1) - (E*(miu/2 - 1/2)*((x - 1)*(y - 1) ...
    + x*y + x*(y - 1) + y*(x - 1) + 2*y*(y - 1)))/(miu^2 - 1) ... 
+ (E*miu*((x - 1)*(y - 1) + x*y + x*(y - 1) + y*(x - 1)))/(miu^2 - 1);

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% mesh generation
n_en   = 4;               % number of nodes in an element
n_el_x = 60;              % number of elements in x-dir
n_el_y = 60;              % number of elements in y-dir
n_el   = n_el_x * n_el_y; % total number of elements

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
    IEN(ee, 1) = (ey-1) * n_np_x + ex;
    IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
    IEN(ee, 3) =  ey    * n_np_x + ex + 1;
    IEN(ee, 4) =  ey    * n_np_x + ex;
  end
end

% ID array
ID = zeros(n_np, 2);
counter = 0;
for ny = 2 : n_np_y - 1
    for nx = 2 : n_np_x - 1
        index = (ny-1)*n_np_x + nx;
        counter = counter + 1;
        ID(index, 1) = counter;
        counter = counter + 1;
        ID(index, 2) = counter;
    end
end

n_eq = counter;

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
for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );
  
  k_ele = zeros(n_en * 2, n_en * 2); % element stiffness matrix
  f_ele = zeros(n_en * 2, 1);    % element load vector
  
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

      f_ele(2*(aa-1)+1) = f_ele(2*(aa-1)+1) + weight(ll) * detJ * f_1(x_l, y_l) * Na;
      f_ele(2*(aa-1)+2) = f_ele(2*(aa-1)+2) + weight(ll) * detJ * f_2(x_l, y_l) * Na;
      
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

% solve the stiffness matrix
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


close all;
[X, Y] = meshgrid(0 : hx : 1, 0 : hy : 1);
Z = reshape(disp(:,2), n_np_x, n_np_y)';
surf(X, Y, Z);

shading interp;

% EOF



% save the solution vector and number of elements to disp with name
% HEAT.mat
save("HEAT", "disp", "n_el_x", "n_el_y");

% EOF

function ei = unit_vector(ii)
    if ii == 1
        ei = [1, 0]';
    else
        ei = [0, 1]';
    end
end