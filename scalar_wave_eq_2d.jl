function scalar_wave_eq_2d_solver(vel, f, source_position, Nx, Ny, h, dt, Nt)
    pml_len = 20;
    pml_coef = 100;
    Nx_pml = Nx + 2*pml_len;
    Ny_pml = Ny + 2*pml_len;
    pml_value = range(0, stop=pml_coef, length=pml_len);
    source_position_pml = source_position .+ pml_len;
    
    u0 = zeros(Nx_pml, Ny_pml);
    u1 = zeros(Nx_pml, Ny_pml);
    u2 = zeros(Nx_pml, Ny_pml);
    wavefield = zeros(Nx, Ny, Nt);
    
    sigma_x, sigma_y = make_pml_coef(Nx, Ny, Nx_pml, Ny_pml, pml_len, pml_value);

    vel_ex = make_vel_ex(vel, Nx_pml, Ny_pml, pml_len);

    # Initialize
    vx1 = zeros(Nx_pml, Ny_pml);
    vx2 = zeros(Nx_pml, Ny_pml);
    vy1 = zeros(Nx_pml, Ny_pml);
    vy2 = zeros(Nx_pml, Ny_pml);
    u0 = zeros(Nx_pml, Ny_pml);
    u1 = zeros(Nx_pml, Ny_pml);
    u2 = zeros(Nx_pml, Ny_pml);
    wavefield = zeros(Nx, Ny, Nt);
    # Coef matrix
    A = ones(Nx_pml, Ny_pml) ./ vel_ex.^2;
    B = (sigma_x + sigma_y) ./ vel_ex.^2;
    C = ones(Nx_pml, Ny_pml) ./ vel_ex.^2;
    # Main loop
    for iter_t = 1:Nt
    
        # original equation
        coef_1 = 2*u1[2:end-1,2:end-1] - u0[2:end-1,2:end-1]  - (dt^2 .* B[2:end-1,2:end-1])./(A[2:end-1,2:end-1]*dt).*(u1[2:end-1,2:end-1]-u0[2:end-1,2:end-1]) - dt^2 .* C[2:end-1,2:end-1]./A[2:end-1,2:end-1].*u1[2:end-1,2:end-1];
    
        coef_2 = dt^2 ./ (A[2:end-1,2:end-1].*(2h)).*(vx1[3:end,2:end-1] - vx1[1:end-2,2:end-1] + vy1[2:end-1,3:end] - vy1[2:end-1,1:end-2]);
    
        coef_3 = dt^2 ./ (A[2:end-1,2:end-1]*h^2).*(u1[3:end,2:end-1] - 2*u1[2:end-1,2:end-1] + u1[1:end-2,2:end-1] + u1[2:end-1,3:end] - 2*u1[2:end-1,2:end-1] + u1[2:end-1,1:end-2]);
    
        u2[2:end-1,2:end-1] = coef_1 + coef_2 + coef_3;
        u2[source_position_pml[1], source_position_pml[2]] += f[iter_t];
    
        # auxiliary equation
        vx2[2:end-1,2:end-1] = vx1[2:end-1,2:end-1] - dt.*sigma_x[2:end-1,2:end-1].*vx1[2:end-1,2:end-1] - dt/(2h).*(sigma_x[2:end-1,2:end-1]-sigma_y[2:end-1,2:end-1]).*(u1[3:end,2:end-1]-u1[1:end-2,2:end-1]);
    
        vy2[2:end-1,2:end-1] = vy1[2:end-1,2:end-1] - dt.*sigma_y[2:end-1,2:end-1].*vy1[2:end-1,2:end-1] - dt/(2h).*(sigma_y[2:end-1,2:end-1]-sigma_x[2:end-1,2:end-1]).*(u1[2:end-1,3:end]-u1[2:end-1,1:end-2]);
    
        # time update
        vx1[:] = vx2; vy1[:] = vy2;
        u0[:] = u1; u1[:] = u2;
    
        # record time domain wavefield
        wavefield[:,:,iter_t] = u2[pml_len+1:pml_len+Nx, pml_len+1:pml_len+Ny];
    end
    return wavefield
end

function make_pml_coef(Nx, Ny, Nx_pml, Ny_pml, pml_len, pml_value)
    sigma_x = zeros(Nx_pml, Ny_pml);
    for i = 1:pml_len
        sigma_x[pml_len+1-i,:] .= pml_value[i];
        sigma_x[pml_len+Nx+i,:] .= pml_value[i];
    end
    sigma_y = zeros(Nx_pml, Ny_pml);
    for i = 1:pml_len
        sigma_y[:,pml_len+1-i] .= pml_value[i];
        sigma_y[:,pml_len+Ny+i] .= pml_value[i];
    end
    return sigma_x, sigma_y
end

function make_vel_ex(vel, Nx_pml, Ny_pml, pml_len)
    vel_ex = zeros(Nx_pml, Ny_pml);
    vel_ex[pml_len+1:end-pml_len, pml_len+1:end-pml_len] .= vel;
    for i = 1:pml_len
        vel_ex[i,:] = vel_ex[pml_len+1,:];
        vel_ex[end-i+1,:] = vel_ex[end-pml_len,:];
        vel_ex[:,i] = vel_ex[:,pml_len+1];
        vel_ex[:,end-i+1] = vel_ex[:,end-pml_len];
    end
    return vel_ex
end

# Ricker function for source
function source_ricker(center_fre, center_time, t)
    x = (1 .- 2*pi^2 .* center_fre^2 .* (t.-center_time).^2) .* exp.(-pi^2*center_fre^2 .* (t .- center_time).^2);
end;