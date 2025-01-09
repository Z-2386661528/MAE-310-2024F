clear all; clc;

run('quarter.m')
transIntoXY;

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

n_g = 5;
[seit, seit_weight] = Gauss(n_g , -1 ,1);


% mesh generation
n_en   = 4;               % number of nodes in an element
n_el   = msh.nbNod;       % total number of elements

% generate the nodal coordinates
xcoor = msh.POS(:,1);
ycoor = msh.POS(:,2);

% IEN array
IEN = msh.QUADS(: , 1:4);

n_lin = size(msh.LINES);
n_lin = n_lin(1,1);
lines_nodel = msh.LINES(:,1:2);
NV = zeros(n_lin,2);

% find the up and right lines
xcounter = 0;nodel_up = zeros(7,3);
ycounter = 0;nodel_right = zeros(7,3);
for ii = 1 : n_lin
    if msh.LINES(ii,3) == 9
        xcounter = xcounter + 1;
        nodel_up(xcounter,:) = msh.LINES(ii,:);
    end
    if msh.LINES(ii,3) == 8
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











