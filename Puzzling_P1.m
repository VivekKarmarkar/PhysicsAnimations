N = 6; % Number of particles
t = 0;
tf = 10;
dt = 0.001;
v_mag = 0.2; % Speed of particles
alpha = 2*pi/N;
t_j = 1;
frames_max = tf/dt + 1;

draw_ngon = true; % switch to draw regular ngon at regular intervals
frames_ngon_gap = 2500;
frames_ngon_max = 3;

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
    particles(j).r(1,:) = [cos(j*alpha) sin(j*alpha) 0];
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
        if t_j < frames_max
            particles(j).r(t_j + 1,:) = particles(j).r(t_j,:) + particles(j).v(t_j,:)*dt;
            particles(j).time(t_j + 1,:) = t + dt;
        end
    end
t = t + dt;
t_j = t_j + 1;
end

figure
xlim([-1.5,1.5])
ylim([-1.5,1.5])
for k=1:N
    r_k = particles(k).r;
    plot(r_k(:,1), r_k(:,2), 'r')
    hold on
    if k<N
        r_k_next = particles(k+1).r;
    else
        r_k_next = particles(1).r;
    end
    line([r_k(1,1) r_k_next(1,1)],[r_k(1,2) r_k_next(1,2)],'color','b')
    if draw_ngon
        for m=1:frames_ngon_max-1
            hold on
            line([r_k(1 + frames_ngon_gap*m,1) r_k_next(1 + frames_ngon_gap*m,1)],[r_k(1 + frames_ngon_gap*m,2) r_k_next(1 + frames_ngon_gap*m,2)],'color','b')
        end
    end
end
