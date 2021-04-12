% Illustrate animation and movies in Matlab
clear
clc
close all

% Step 1: Load data
load('particles.mat')
N = 6; % Need to grab N from particles or make all of this into one code
t = particles(1).time;
frame_count = 0;

% Step 2: Draw/Render the scenario
figh = figure;

for k=1:100:length(t)
    % Wipe the slate clean so that we are plotting with a blank figure
    clf
    
    for j=1:N
        % Generate the data for current particle
        t = particles(j).time;
        x = particles(j).r(:,1);
        y = particles(j).r(:,2);
        % Extract the data at the current time
        t_k = t(k);
        x_k = x(k);
        y_k = y(k);

        % Plot the current location of the particle
        plot(x_k, y_k, 'go', 'LineWidth', 3, 'MarkerSize', 15)

        % Plot the entire trajectory
        hold on
        plot(x, y, 'b-', 'LineWidth', 2)
    end
    
    % Decorate the plot
    grid on
    xlabel('x')
    ylabel('y')
    
    % Save the image information regarding each frame
    frame_count = frame_count + 1;
    movie_vector(frame_count) = getframe(gcf);
end

% Save the movie
my_writer = VideoWriter('Puzzling200P1Vid', 'MPEG-4');
my_writer.FrameRate = 20;

% Open the VideoWriter object, write the movie and close the file
open(my_writer);
writeVideo(my_writer, movie_vector);
close(my_writer);