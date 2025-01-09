clear all; clc;

run('quarter.m')
transIntoXY;

f_1 = @(x,y) 0;f_2 = @(x,y) 0;

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
n_eq = 224;
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

% generate the nodal coordinates
xcoor = msh.POS(:,1);
ycoor = msh.POS(:,2);

% IEN array
IEN = msh.QUADS(: , 1:4);

% 右边和上边条件设置
Dirichlet = 0;
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
    if ID(ii,1) ~= 0
        ID(ii,1) = counter;
        counter = counter + 1;
    elseif ID(ii,2) ~= 0
        ID(ii,2) = counter;
        counter = counter + 1;
    end
end

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
xcounter = 0;nodel_up = zeros(7,3);
ycounter = 0;nodel_right = zeros(7,3);
for ii = 1 : n_lin
    if msh.LINES(ii,3) == 9
        xcounter = xcounter + 1;
        nodel_up(xcounter,:) = msh.LINES(ii,:);
    elseif msh.LINES(ii,3) == 8
        ycounter = ycounter + 1;
        nodel_right(ycounter,:) = msh.LINES(ii,:);
    end
end
% for gauss point on lines and get the sigma
xpos = zeros(7,2);ypos = zeros(7,2);
gaussPointX = zeros(7,n_g);gaussPointY = zeros(7,n_g);
upper_sigma = cell(7,n_g);right_sigma = cell(7,n_g);

for ii = 1 : 7
    xpos(ii ,:) = msh.POS(nodel_up(ii,1:2));
    right_nodel_list = msh.POS(nodel_right(ii,1:2),:);
    ypos(ii ,:) = right_nodel_list(:,2);

    [upper_nx(ii,1) , upper_nx(ii,2)] = normal_vector(xpos(ii ,1),xpos(ii ,2),1,1); % find the normal vector
    [upper_ny(ii,1) , upper_ny(ii,2)] = normal_vector(1,1,ypos(ii ,1),ypos(ii ,2));

    for jj = 1 : n_g
        dertx = xpos(ii ,1) - xpos(ii ,2);
        gaussPointX(ii,jj) = xpos(ii ,1) - dertx *(1 - seit(jj))/2;
        derty = ypos(ii ,1) - ypos(ii ,2);
        gaussPointY(ii,jj) = ypos(ii ,1) - derty *(1 - seit(jj))/2;

        upper_sigma{ii , jj} = zeros(2,2);
        upper_sigma{ii , jj}(1,1) = subs(sigma_xx, [x, y], [gaussPointX(ii , jj) + 1, 2]);
        upper_sigma{ii , jj}(2,2) = subs(sigma_yy, [x, y], [gaussPointX(ii , jj) + 1, 2]);
        upper_sigma{ii , jj}(1,2) = subs(tao_xy, [x, y], [gaussPointX(ii , jj) + 1, 2]);
        upper_sigma{ii , jj}(2,1) = upper_sigma{ii , jj}(1,2);

        right_sigma{ii , jj} = zeros(2,2);
        right_sigma{ii , jj}(1,1) = subs(sigma_xx, [x, y], [2,gaussPointY(ii , jj) + 1]);
        right_sigma{ii , jj}(2,2) = subs(sigma_yy, [x, y], [2,gaussPointY(ii , jj) + 1]);
        right_sigma{ii , jj}(1,2) = subs(tao_xy, [x, y], [2,gaussPointY(ii , jj) + 1]);
        right_sigma{ii , jj}(2,1) = right_sigma{ii , jj}(1,2);
    end
end

upper_h = zeros(n_g * 7,2);         % Calculate hi
right_h = zeros(n_g * 7,2);
encounter = 0;
for ii = 1 : 7
    for jj = 1 : n_g
        encounter = encounter + 1;
        upper_h(encounter,:) = upper_sigma{ii,jj} * [upper_nx(ii,1) , upper_nx(ii,2)]';
        right_h(encounter,:) = right_sigma{ii,jj} * [upper_ny(ii,1) , upper_ny(ii,2)]';
    end
end

% allocate the stiffness matrix and load vector
% K = spalloc(n_eq, n_eq, 9 * n_eq);
K = zeros(n_eq, n_eq);
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










function ei = unit_vector(ii)
    if ii == 1
        ei = [1, 0]';
    else
        ei = [0, 1]';
    end
end