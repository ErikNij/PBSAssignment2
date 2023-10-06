% Read data from CSV file
data = readtable('Save_Energy.csv');

% Extract Time and energy values
Time = data.Time;
Epot = data.Epot;
Ekin = data.Ekin;
Etot = data.Etot;

% Plotting the data
figure;
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

% Display the plot
