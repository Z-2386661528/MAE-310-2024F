clear all; clc;

run('quarter.m')

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

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
line_POS = cell(2, 1);
NV = zeros(n_lin,2);

% Normal Vector of lines
for ii = 1 : n_lin
    for jj = 1 : 2
        for kk = 1 : 2
            line_POS{kk}(ii,jj) = msh.POS(lines_nodel(ii,jj),kk);
        end
    end
end

for ii = 1 : n_lin
    norm = 0.0;
    NV(ii,1) = line_POS{2}(ii,1) - line_POS{2}(ii,2);
    NV(ii,2) = - line_POS{1}(ii,1) + line_POS{1}(ii,2);
    norm = (NV(ii,1)^2 + NV(ii,2)^2)^0.5;
    NV(ii,1) = NV(ii,1) / norm;
    NV(ii,2) = NV(ii,2) / norm;    
end
