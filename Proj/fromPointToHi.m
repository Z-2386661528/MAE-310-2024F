function [hx, hy] = fromPointToHi(x , y , nx , ny)

T = 1E4;
R = 0.5;

r = sqrt(x^2+y^2);
seit = atan2(y,x);

sigma_rr = T / 2 * (1 - (R/r) ^ 2) + T / 2 * (1 - 4 * (R/r)^2 + 3 *(R/r)^4)*cos(2*seit);
sigma_tt = T / 2 * (1 + (R/r) ^ 2) - T / 2 * (1 + 3*(R/r)^4)*cos(2*seit);
tao_rt   = - T / 2 * (1 + 2 * (R/r)^2 - 3 *(R/r)^4)*sin(2*seit);

sigma_xx = sigma_rr * cos(-seit)^2 + sigma_tt * sin(-seit)^2 + sin(-2*seit) * tao_rt;
sigma_yy = sigma_rr * sin(-seit)^2 + sigma_tt * cos(-seit)^2 - sin(-2*seit) * tao_rt;
tao_xy   = - sigma_rr * sin(-2*seit) / 2 + sigma_tt * sin(-2*seit) / 2 + tao_rt*cos(-2 * seit);



sigma_xxx = @(x,y) sigma_xx;
sigma_yyx = @(x,y) sigma_yy;
tao_xyx = @(x,y) tao_xy;

sigma = zeros(2,2);
sigma(1,1) = sigma_xxx(x+1,y+1);
sigma(2,2) = sigma_yyx(x+1,y+1);
sigma(2,1) = tao_xyx(x+1,y+1);
sigma(1,2) = sigma(2,1);

h = sigma * [nx , ny]';
hx=h(1,1);
hy=h(2,1);