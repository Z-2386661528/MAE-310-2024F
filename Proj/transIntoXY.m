syms x y
T = 1;
R = 0.5;

r = sqrt(x^2+y^2);
seit = atan2(y,x);

sigma_rr = T / 2 * (1 - (R/r) ^ 2) + T / 2 * (1 - 4 * (R/r)^2 + 3 *(R/r)^4)*cos(2*seit);
sigma_tt = T / 2 * (1 + (R/r) ^ 2) - T / 2 * (1 + 3*(R/r)^4)*cos(2*seit);
tao_rt   = - T / 2 * (1 + 2 * (R/r)^2 - 3 *(R/r)^4)*sin(2*seit);

sigma_xx = sigma_rr * cos(-seit)^2 + sigma_tt * sin(-seit)^2 + sin(-2*seit) * tao_rt;
sigma_yy = sigma_rr * sin(-seit)^2 + sigma_tt * cos(-seit)^2 - sin(-2*seit) * tao_rt;
tao_xy   = - sigma_rr * sin(-2*seit) / 2 + sigma_tt * sin(-2*seit) / 2 + tao_rt*cos(-2 * seit);

sigma_xx = simplify(sigma_xx);
sigma_yy = simplify(sigma_yy);
tao_xy   = simplify(tao_xy);
