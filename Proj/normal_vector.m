function [nx, ny] = normal_vector(x_1,x_2,y_1,y_2)

dertx = x_2 - x_1;
derty = - y_2 + y_1;
norm = (dertx^2 + derty^2)^0.5;
nx = derty / norm;
ny = dertx / norm; 

