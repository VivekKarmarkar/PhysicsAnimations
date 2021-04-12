clear
clc
close all

N = 6; % Number of particles
t = 0;
r_ini = 1;
dt = 0.001;
v_mag = 0.2; % Speed of particles
alpha = 2*pi/N;
d_ini = 2*r_ini*sin(alpha/2);
tf = d_ini/(v_mag*(1-cos(alpha))); % calculate relative velocity by considering component velocities of both particles
t_j = 1;
frames_max = int16(tf/dt + 1);

asymmetry = true; % Breaking the symmetry of the problem by giving Nth particle twice the speed

particles = struct;
particles(N).time = zeros(frames_max,1);
particles(N).r = zeros(frames_max,3);
particles(N).r_rel = zeros(frames_max,3);
particles(N).v = v_mag*ones(frames_max,3);
particles(N).n_hat = ones(frames_max,3);
for j=1:N
    particles(j).r = zeros(frames_max,3);
    particles(j).r_rel = zeros(frames_max,3);
    particles(j).v = v_mag*ones(frames_max,3);
    particles(j).n_hat = ones(frames_max,3);
    particles(j).r(1,:) = [r_ini*cos(j*alpha) r_ini*sin(j*alpha) 0];
    particles(j).time(1) = 0;
end

while t<tf
    for j=1:N
        if j<N
            particles(j).r_rel(t_j,:) = particles(j+1).r(t_j,:) - particles(j).r(t_j,:);
        else
            particles(j).r_rel(t_j,:) = particles(1).r(t_j,:) - particles(j).r(t_j,:);
        end
        particles(j).n_hat(t_j,:) = particles(j).r_rel(t_j,:) ./ norm(particles(j).r_rel(t_j,:));
        particles(j).v(t_j,:) = v_mag*particles(j).n_hat(t_j,:);
        if asymmetry
            particles(N).v(t_j,:) = 2*v_mag*particles(N).n_hat(t_j,:);
        end
        if t_j < frames_max
            particles(j).r(t_j + 1,:) = particles(j).r(t_j,:) + particles(j).v(t_j,:)*dt;
            particles(j).time(t_j + 1,:) = t + dt;
        end
    end
t = t + dt;
t_j = t_j + 1;
end

t_new = particles(1).time;
frame_count = 0;

disp("WARNING: DO NOT CLOSE THE ANIMATION while MATLAB reads BUSY in left bottom corner")

figh = figure;

for k=1:100:length(t_new)
    clf
    
    for j=1:N
        T = particles(j).time;
        x = particles(j).r(:,1);
        y = particles(j).r(:,2);
        if j<N
            x_next = particles(j+1).r(:,1);
            y_next = particles(j+1).r(:,2);
        else
            x_next = particles(1).r(:,1);
            y_next = particles(1).r(:,2);
        end
        
        t_k = T(k);
        x_k = x(k);
        y_k = y(k);
        
        x_k_next = x_next(k);
        y_k_next = y_next(k);

        plot(x_k, y_k, 'go', 'LineWidth', 3, 'MarkerSize', 15)
        hold on
        line([x_k x_k_next],[y_k y_k_next], 'color', 'r', 'LineWidth', 2);

        hold on
        plot(x, y, 'b-', 'LineWidth', 2)
    end
    
    grid on
    xlabel('x')
    ylabel('y')
    
    frame_count = frame_count + 1;
    movie_vector(frame_count) = getframe(gcf);
end

my_writer = VideoWriter('Puzzling200P1Vid', 'MPEG-4');
my_writer.FrameRate = 20;

open(my_writer);
writeVideo(my_writer, movie_vector);
close(my_writer);
