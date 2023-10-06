% Read the data from the CSV file
data = csvread('Velocity_Distribution.csv', 1, 0);

% Extract velocity and probability values
velocity = data(:, 1);
probability = data(:, 2);

% Plot the probability density function
figure;
plot(velocity, probability, 'LineWidth', 1.5);
xlabel('Velocity');
ylabel('Probability Density');
title('Velocity Probability Density Function');
grid on;


