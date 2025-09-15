%% Problem 5

clear; clc; close all;

[t_vec, av_pos_inert, av_att, tar_pos_inert, tar_att] = LoadASPENData('3801_Sec001_Test1.csv');

n = length(t_vec);

av_att313 = zeros(3,n);
tar_att313 = zeros(3,n);

for k = 1:n
    C_EB_av = RotationMatrix321(av_att(:,k));   
    av_att313(:,k) = EulerAngles313(C_EB_av);  
    
    C_EB_tar = RotationMatrix321(tar_att(:,k));
    tar_att313(:,k) = EulerAngles313(C_EB_tar);
end

av_att313_deg = rad2deg(av_att313);
tar_att313_deg = rad2deg(tar_att313);

figure;

subplot(3,1,1);
plot(t_vec, av_att313_deg(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(t_vec, tar_att313_deg(1,:), 'r--', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('\alpha [deg]');
title('Problem 5: 3-1-3 Euler Angle \alpha');
legend('Aerospace Vehicle','Target');
grid on;

subplot(3,1,2);
plot(t_vec, av_att313_deg(2,:), 'b', 'LineWidth', 1.5); hold on;
plot(t_vec, tar_att313_deg(2,:), 'r--', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('\beta [deg]');
title('Problem 5: 3-1-3 Euler Angle \beta');
grid on;

subplot(3,1,3);
plot(t_vec, av_att313_deg(3,:), 'b', 'LineWidth', 1.5); hold on;
plot(t_vec, tar_att313_deg(3,:), 'r--', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('\gamma [deg]');
title('Problem 5: 3-1-3 Euler Angle \gamma');
grid on;

sgtitle('Problem 5: 3-1-3 Euler Angles vs Time (Vehicle: blue, Target: red)');


%% Problem 6
rel_pos_E = tar_pos_inert - av_pos_inert; 
figure;

subplot(3,1,1);
plot(t_vec, rel_pos_E(1,:), 'b', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('x_E [m]');
title('Problem 6: Relative Position in Frame E (x-component)');
grid on;

subplot(3,1,2);
plot(t_vec, rel_pos_E(2,:), 'r', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('y_E [m]');
title('Problem 6: Relative Position in Frame E (y-component)');
grid on;

subplot(3,1,3);
plot(t_vec, rel_pos_E(3,:), 'g', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('z_E [m]');
title('Problem 6: Relative Position in Frame E (z-component)');
grid on;

sgtitle('Problem 6: Relative Position of Target w.r.t Aerospace Vehicle in Frame E');


%% Problem 7
rel_pos_B = zeros(3,n);

for k = 1:n
   
    C_EB = RotationMatrix321(av_att(:,k));   
    C_BE = C_EB.';                           
    
    rel_pos_B(:,k) = C_BE * rel_pos_E(:,k);
end

figure;

subplot(3,1,1);
plot(t_vec, rel_pos_B(1,:), 'b', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('x_B [m]');
title('Problem 7: Relative Position in Frame B (x-component)');
grid on;

subplot(3,1,2);
plot(t_vec, rel_pos_B(2,:), 'r', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('y_B [m]');
title('Problem 7: Relative Position in Frame B (y-component)');
grid on;

subplot(3,1,3);
plot(t_vec, rel_pos_B(3,:), 'g', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('z_B [m]');
title('Problem 7: Relative Position in Frame B (z-component)');
grid on;

sgtitle('Problem 7: Relative Position of Target w.r.t Aerospace Vehicle in Frame B');
