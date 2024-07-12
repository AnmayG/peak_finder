function peak_freqs = params_to_peaks(params)
    phi = params.phi * pi / 180;
    theta = params.theta * pi / 180;
    phi_E = params.phi_E * pi / 180;
    theta_E = params.theta_E * pi / 180;
    sxxpyy = params.sxxpyy;
    sxxnyy = params.sxxnyy;
    szz = params.szz;
    sxz = params.sxz;
    syz = params.syz;
    sxy = params.sxy;

    stress_tensor = [[(sxxpyy + sxxnyy) / 2, sxy, sxz]; ...
                     [sxy, (sxxpyy - sxxnyy) / 2, syz]; ...
                     [sxz, syz, szz]];
    b_direction = [cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)];
    e_direction = [cos(phi_E) * sin(theta_E), sin(phi_E) * sin(theta_E), cos(theta_E)];

    % Calculate axis directions
    zz = zeros(4, 3);
    if params.orientation
        % Using 100 orientation where all groups can be clearly separated
        zz(1, :) = [1, 1, 1] / sqrt(3);
        zz(2, :) = [-1, -1, 1] / sqrt(3);
        zz(3, :) = [-1, 1, -1] / sqrt(3);
        zz(4, :) = [1, -1, -1] / sqrt(3);
    else
        % Using 111 orientation where there's one clear axis to align with
        zz(1, :) = [0, 0, 1];
        zz(2, :) = [2 * sqrt(2) / 3, 0, -1/3];
        zz(3, :) = [-sqrt(2) / 3, sqrt(6)/3, -1/3];
        zz(4, :) = [-sqrt(2) / 3, -sqrt(6)/3, -1/3];
    end

    xx = zeros(4, 3);
    for i = 1:4
        zz_ind = mod(i, 4) + 1;
        projection = zz(zz_ind, :) - dot(zz(i, :), zz(zz_ind, :)) .* zz(i, :);
        xx(i, :) = projection / norm(projection);
    end

    yy = zeros(4, 3);
    for i = 1:4
        yy(i, :) = cross(zz(i, :), xx(i, :)) / norm(cross(zz(i, :), xx(i, :)));
    end

    % Make a list of rotation matrices to rotate between the different
    % frames of the stress tensor
    R = zeros(3, 3, 4);
    for i=1:4
        mm = [xx(i, :); yy(i, :); zz(i, :)];
        R(:, :, i) = mm';
    end

    alpha = zeros(4);
    beta = zeros(4);
    alpha_E = zeros(4);
    beta_E = zeros(4);
    % Spin operators
    Sx = [[0, 1, 0]; [1, 0, 1]; [0, 1, 0]]/sqrt(2);
    Sy = [[0, 1, 0]; [-1, 0, 1]; [0, -1, 0]]/(sqrt(2)*1i);
    Sz = [[1, 0, 0]; [0, 0, 0]; [0, 0, -1]];
    
    for i=1:4
        % angle between NV axis and B field
        b_angle = dot(b_direction, zz(i, :));
        alpha(i) = acos(b_angle);
        % angle between x-axis and non-NV-aligned components of B field
        beta(i) = acos(dot((b_direction - b_angle) / norm(b_direction - b_angle), xx(i, :)));
        e_angle = dot(e_direction, zz(i, :));
        alpha_E(i) = acos(e_angle);
        beta_E(i) = acos(dot((e_direction - e_angle) / norm(e_direction - e_angle), xx(i, :)));
    end

    % Stress coupling constants
    a1 = 4.86*1e6;
    a2 = -3.7*1e6;        
    b = -2.3*1e6;
    c = 3.5*1e6;

    % Magnetic and electric field constants
    gamma = 28025 * 100;
    parE = 0;
    perpE = 17;
    % Solve for the Hamiltonian to get whatever groups we need
    peak_freqs = zeros(8);
    for i=1:4 % Iterate through all of the axes and compute the Hamiltonian
        rot = squeeze(R(:, :, i));
        S = rot * stress_tensor * rot';
        xx = S(1, 1);
        xy = S(1, 2);
        xz = S(1, 3);
        yy = S(2, 2);
        yz = S(2, 3);
        zz = S(3, 3);

        Mz = (a1-a2)*(xx+yy)+(a1+2*a2)*zz;  % a1*(xx+yy+zz)+2*a2*(yz+xz+xy)
        Mx = -(b+c)*(yy-xx)+(sqrt(2)*b-sqrt(2)/2*c)*2*xz;
        My = -(b+c)*2*xy+(sqrt(2)*b-sqrt(2)/2*c)*2*yz;
        
        zfs_term = params.ZFS * dot(Sz, Sz);
        magnetic_term = gamma * params.BB * ... 
            (cos(alpha(i)) * Sz + ...
            sin(alpha(i)) * cos(beta(i)) * Sx + ...
            sin(alpha(i)) * sin(beta(i)) * Sy);
        electric_term = params.EE * (parE * sin(theta_E) * dot(Sz, Sz) + ...
                             perpE * cos(theta_E) * cos(phi_E) * (dot(Sy, Sy)-dot(Sx, Sx)) + ...
                             perpE * cos(theta_E) * sin(phi_E) * (dot(Sx, Sy) + dot(Sy, Sx)));
        spin_term = Mz * dot(Sz, Sz) + Mx * (dot(Sy, Sy)-dot(Sx, Sx)) + My * (dot(Sx, Sy) + dot(Sy, Sx));
        HH = zfs_term + magnetic_term + electric_term + spin_term;
        eigval = sort(real(eig(HH)));
        fl = eigval(2) - eigval(1);
        fh = eigval(3) - eigval(1);
        peak_freqs(2 * i - 1) = fl;
        peak_freqs(2 * i) = fh;
    end
    disp(peak_freqs)
end