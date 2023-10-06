%% energy
% Read data from CSV file
data1 = readtable('Save_Energy.csv');

% Extract Time and energy values
Time = data1.Time;
Epot = data1.Epot;
Ekin = data1.Ekin;
Etot = data1.Etot;

% Plotting the data
figure(1);
plot(Time, Epot, 'b-', 'LineWidth', 2, 'DisplayName', 'Epot');
hold on;
plot(Time, Ekin, 'r-', 'LineWidth', 2, 'DisplayName', 'Ekin');
plot(Time, Etot, 'g-', 'LineWidth', 2, 'DisplayName', 'Etot');
hold off;

% Add labels, title, and legend
xlabel('Time [seconds]');
ylabel('Energy [uA^2fs^{-2}]');
title('Energy vs Time');
legend('Location', 'Best');
grid on;


%% pdf
% Read the data from the CSV file
data2 = csvread('Velocity_Distribution.csv', 1, 0);

% Extract velocity and probability values
velocity = data2(:, 1);
probability = data2(:, 2);

% Plot the probability density function
figure(2);
plot(velocity, probability, 'LineWidth', 1.5);
xlabel('Velocity');
ylabel('Probability Density');
title('Velocity Probability Density Function');
grid on;
