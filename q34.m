%% Problem 3
% Figure with the 3d positions of the aerospace vehicle and the target in
% Frame N

figure('units','inches','Position', [2 2 7 5]) 
hold on
plot3(pos_av_aspen(1,:), pos_av_aspen(2,:), pos_av_aspen(3,:), 'b') % Aerospace vehicle position in Frame N
plot3(pos_tar_aspen(1,:), pos_tar_aspen(2,:), pos_tar_aspen(3,:), 'r--') % Target position in Frame N
legend('Aerospace Vehicle Position', 'Target Position')
title('Problem 3, Position of the Aerospace Vehicle Compared to the Target in the N Frame')
xlabel('X Position (m)')
ylabel('Y Position (m)')
zlabel('Z Position (m)')
view(23.8956521270877, 19.4042801275688)
print('Q3', '-dpng');

%% Problem 4
% Figure for the subplot of the changes in x,y,z position of objects in the
% E frame
figure(2)
sgtitle('Problem 4p1, Position Vector Components in Frame E vs Time')
subplot(3, 1, 1)
hold on
plot(t_vec, av_pos_inert(1,:))
plot(t_vec, tar_pos_inert(1,:), 'r--')
title('X Position in Frame E vs Time')
ylabel('X Position (m)')
xlabel('Time (s)')

subplot(3, 1, 2)
hold on
plot(t_vec, av_pos_inert(2,:), 'b')
plot(t_vec, tar_pos_inert(2,:), 'r--')
title('Y Position in Frame E vs Time')
ylabel('Y Position (m)')
xlabel('Time (s)')

subplot(3, 1, 3)
hold on
plot(t_vec, av_pos_inert(3,:), 'b')
plot(t_vec, tar_pos_inert(3,:), 'r--')
title('Z Position in Frame E vs Time')
ylabel('Z Position (m)')
xlabel('Time (s)')

fig = gcf;
fig.Position(4) = fig.Position(4) + 80;
lgnd = legend('Aerospace Vehicle', 'Target', 'show');
lgnd.Position(1) = 0.05;
lgnd.Position(2) = 0.02;
print('Q4p1', '-dpng');

figure(3)
sgtitle('Problem 4p2, Euler Angle Components in Frame E v Time')
subplot(3, 1, 1) % Roll
hold on
plot(t_vec, av_att(1, :)*180/pi, 'b')
plot(t_vec, tar_att(1,:)*180/pi, 'r--')
title('Roll in Frame E vs Time')
ylabel('Roll Angle (deg)')
xlabel('Time (s)')

subplot(3, 1, 2) % Pitch
hold on
plot(t_vec, av_att(2, :)*180/pi, 'b')
plot(t_vec, tar_att(2,:)*180/pi, 'r--')
title('Pitch in Frame E vs Time')
ylabel('Pitch Angle (deg)')
xlabel('Time (s)')

subplot(3, 1, 3) % Yaw
hold on
plot(t_vec, av_att(3, :)*180/pi, 'b')
plot(t_vec, tar_att(3,:)*180/pi, 'r--')
title('Yaw in Frame E vs Time')
ylabel('Yaw Angle (deg)')
xlabel('Time (s)')

fig = gcf;
fig.Position(4) = fig.Position(4) + 80;
lgnd = legend('Aerospace Vehicle', 'Target', 'show');
lgnd.Position(1) = 0.05;
lgnd.Position(2) = 0.02;
print('Q4p2', '-dpng');
