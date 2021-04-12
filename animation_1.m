% Illustrate animation and movies in Matlab
clear
clc
close all

% Step 1: Generate data
t = linspace(0, 2*pi, 100);
x = 5*cos(t);
y = 2*sin(t);
z = t;

% Step 2: Draw/Render the scenario
figh = figure;

for k=1:length(t)
    % Wipe the slate clean so that we are plotting with a blank figure
    clf
    
    % Extract the data at the current time
    t_k = t(k);
    x_k = x(k);
    y_k = y(k);
    z_k = z(k);
    
    % Plot the current location of the particle
    plot3(x_k, y_k, z_k, 'go', 'LineWidth', 3, 'MarkerSize', 15)
    
    % Plot the entire trajectory
    hold on
    plot3(x, y, z, 'b-', 'LineWidth', 2)
    
    % Decorate the plot
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title(['Particle at t = ', num2str(t_k), ' seconds'])
    view ([30 + 20*t_k 35]) % Rotate the field of view
    
    % Save the image information regarding each frame
    movie_vector(k) = getframe(figh, [10 10 520 400]); % Use standardized dimensions
end

% Save the movie
my_writer = VideoWriter('curve', 'MPEG-4');
my_writer.FrameRate = 20;

% Open the VideoWriter object, write the movie and close the file
open(my_writer);
writeVideo(my_writer, movie_vector);
close(my_writer);