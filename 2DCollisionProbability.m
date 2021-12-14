rho_ab_1 = 0.2;
rho_bc_1 = -0.5;
rho_ac_1 = -0.8;
sigma_a_1 = 5;
sigma_b_1 = 4;
sigma_c_1 = 10;
C1 = [sigma_a_1^2 rho_ab_1*sigma_a_1*sigma_b_1 rho_ac_1*sigma_a_1*sigma_c_1; rho_ab_1*sigma_a_1*sigma_b_1 sigma_b_1^2 rho_bc_1*sigma_b_1*sigma_c_1; rho_ac_1*sigma_a_1*sigma_c_1 rho_bc_1*sigma_b_1*sigma_c_1 sigma_c_1^2];

rho_ab_2 = -0.3;
rho_bc_2 = -0.1;
rho_ac_2 = 0.4;
sigma_a_2 = 1;
sigma_b_2 = 2;
sigma_c_2 = 3;
C2 = [sigma_a_2^2 rho_ab_2*sigma_a_2*sigma_b_2 rho_ac_2*sigma_a_2*sigma_c_2; rho_ab_2*sigma_a_2*sigma_b_2 sigma_b_2^2 rho_bc_2*sigma_b_2*sigma_c_2; rho_ac_2*sigma_a_2*sigma_c_2 rho_bc_2*sigma_b_2*sigma_c_2 sigma_c_2^2];

r1  = [378.39559 4305.721887 5752.767554];
v1  = [2.360800244 5.580331936 -4.322349039];

r2  = [374.5180598 4307.560983 5751.130418];
v2  = [-5.388125081 -3.946827739 3.322820358];

Rel_Conv = 1e-09;
HBR = 0.02;
HBR_type = 'square';

C_total = C1 + C2;
% relative reference frame
r = r1-r2;
v = v1-v2;
h = cross(r,v);
% from definition of frame
y = v/norm(v);
z = h/norm(h);
x = cross(y,z);

transformation = [x;y;z];
trans_cov_total = transformation * C_total * transpose(transformation);
twod_proj = [1 0 0; 0 0 1];
C = twod_proj * trans_cov_total * transpose(twod_proj);

if ~chol(C)==0
    error('please check covariance matrix')
end

% HBR is center of encounter plane
x0 = norm(r);
z0 = 0;
C_inv = inv(C);

% Integrand
r_proj = [x ; z]
Integrand = exp(-1/2*r_proj'*C_inv*r_proj)

Abs_Conv = 1e-13;
    switch HBR_type
        case 'square'
            % define semi-circle w.r.t x as the variable 
            top_half = @(x)(sqrt(HBR^2 - (x-x0).^2) .* (abs(x-x0)<=HBR));
            lower_half = @(x)(-sqrt(HBR^2 - (x-x0).^2) .* (abs(x-x0)<=HBR));
            Pc = 1/(2*pi*sqrt(det(C)))*quad2d(Integrand,x0-HBR,x0+HBR,top_half,lower_half,'RelTol',Rel_Conv,'AbsTol',Abs_Conv);
        case 'circle'
            Pc = 1/(2*pi)*1/sqrt(det(C))*quad2d(Integrand,x0-HBR,x0+HBR,z0-HBR,z0+HBR,'AbsTol',Abs_Conv,'RelTol',Rel_Conv);
    end





