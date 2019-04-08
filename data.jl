Nx = 101; Ny = 101;
dx = 0.005; dy = dx;

Fs = 1000;
dt = 1/Fs;
Nt = 701;
t = range(0, step=dt, length=Nt);

# source_func = sin.(2*pi*10*t);
source_func = source_ricker(10, 0.1, t);
# plot(t, source_func); title("Source Function")

source = zeros(Nt,6);
for i = 1:6
    source[:,i] .= source_func;
end
source[:,6] .=  sin.(2*pi*10*t);

source_position = zeros(Int, 6, 2);
for i = 1:6
    source_position[i,2] = (i-1)*20 + 1
    source_position[i,1] = 3;
end

# Model 1
receiver_position = zeros(Int,101,2);
for i = 1:101
    receiver_position[i,2] = (i-1)*1 + 1
    receiver_position[i,1] = 1;
end

c = ones(Nx, Ny); rho = ones(Nx,Ny);
a0 = 1 ./ (c.^2 .* rho);
b0 = 1 ./ rho;

c[50:end,:] .= 1.2;

a = 1 ./ (c.^2 .* rho);
b = 1 ./ rho;
# using ImageFiltering
# a0 = imfilter(a, Kernel.gaussian(30));
# b0 = imfilter(b, Kernel.gaussian(30));