%% Read OptiTrack File
clear; close all; clc;
filename = 'me133a_6nov_wand_250_270m.csv';
fileraw = csvread(filename,7,0);
time = fileraw(:,2);
rb_q = fileraw(:,3:6)';
rb_xyz = (fileraw(:,7:9))';
nmarkers = 5;
marker = cell(nmarkers,1);
for i=1:nmarkers
   marker{i} =  fileraw(:,7+4*i:9+4*i);
end

%% Rodriguez Screw-Parameters Approach
nmarkers = 3;
ndifferences = length(marker{1})-1;

% Construct Q and P matrices
PQR0 = cell(nmarkers,1);
PQR1 = cell(nmarkers,1);
for i=1:nmarkers
   PQR0{i} = marker{i}(1:end-1,:);
   PQR1{i} = marker{i}(2:end  ,:);
end

% compute vector tan(phi/2)*w
tw = zeros(nmarkers-1,3);
for i = 1:ndifferences
    % define P, Q, R 0 and 1 for easy calculation
    P0 = PQR0{1}(i, :);
    P1 = PQR1{1}(i, :);
    Q0 = PQR0{2}(i, :);
    Q1 = PQR1{2}(i, :);
    R0 = PQR0{3}(i, :);
    R1 = PQR1{3}(i, :);
    
    tw(i,:) = cross((Q1 - Q0) - (R1 - R0), (P1 - P0) - (R1 - R0)) / ...
              dot((Q1 - Q0) - (R1 - R0), (P1 + P0) - (R1 + R0));
end

% take the magnitude to get phi
tanphi2 = vecnorm(tw,2,2);
phi = 2*atan(tanphi2);  % in radians, may have sign error due to atan

% get w (omega) by making tw into unit vector
w = tw ./ tanphi2;

% rhoperp is easy to get
rhoperp = zeros(nmarkers-1,3);
for i = 1:ndifferences
    % define P, Q, R 0 and 1 for easy calculation
    P0 = PQR0{1}(i, :);
    P1 = PQR1{1}(i, :);
    Q0 = PQR0{2}(i, :);
    Q1 = PQR1{2}(i, :);
    R0 = PQR0{3}(i, :);
    R1 = PQR1{3}(i, :);
    
    rhoperp(i,:) = 0.5*(cross(w(i,:), P1 - P0) ./ tanphi2(i) - dot(w(i,:), P1 + P0)*w(i,:) + P0 + P1);
end

% d parallel (dpar) is even easier
dpar = zeros(nmarkers-1,1);
for i = 1:ndifferences
    % just need p0 p1
    P0 = PQR0{1}(i, :);
    P1 = PQR1{1}(i, :);
    
    dpar(i,:) = dot(w(i,:), P1 - P0);
end

% pitch (h) = dpar / phi
h = dpar ./ phi;

%% compute quaternion
q = [cos(phi/2), w(:,1).*sin(phi/2), w(:,2).*sin(phi/2), w(:,3).*sin(phi/2)];
q = q'; % transpose to match rb_q

% rearrange rb_q to match
rb_q = [ rb_q(4,:); rb_q(1,:); rb_q(2,:); rb_q(3,:)];

%% Find displacements between each adjacent pair of frames i (in rb_q)
% for the provided optitrack data
opti_change = zeros(size(q));
for j = 1:size(rb_q,2)-1
    qk = rb_q(:,j+1);
    qi = rb_q(:,j);
    opti_change(:,j) = quatmultiply(qk', quatinv(qi'));
end
opti_change(:,1200) = NaN;

% Now instead of successive frames we compare to 1st frame
% We convert our quaternion displacement back
q2 = zeros(size(rb_q));
q2(:,1) = rb_q(:,1);
for j = 1:size(rb_q,2)-1
    q2(:,j+1) = quatmultiply(q(:,j)', q2(:,j)');
end

% q stores quaternion displacement between successive frames
% q2 stores quaternion displacement relative to 1st frame

%% Compute displacement
% Get the skew matrix corresponding to w (what)
what = zeros(3, 3, ndifferences);
for i = 1:ndifferences
    what(:, :, i) = [      0, -w(i,3),  w(i,2); ...
                      w(i,3),       0, -w(i,1); ...
                     -w(i,2),  w(i,1),      0];
end

% Convert these parameters into homogeneous coordinates
% R_AB
R = zeros(3, 3, ndifferences);
for i = 1:ndifferences
    R(:,:,i) = expm(phi(i)*what(:,:,i)); 
end

% check R
R131 = zeros(3,3,ndifferences);
for i = 1:ndifferences
    a = q(1, i);
    b = q(2, i);
    c = q(3, i);
    d = q(4, i);
    R131(:,:,i) = [[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d + a*c)];...
        [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)];...
        [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]]; 
end

% d_AB
d = zeros(ndifferences, 3);
for i = 1:ndifferences
    d(i,:) = ((eye(3) - R(:,:,i))*rhoperp(i,:)')' + h(i)*phi(i)*w(i,:); 
end

% put all together into g_AB
gab = zeros(4, 4, ndifferences);
for i = 1:ndifferences
    Ri = R(:,:,i);
    di = d(i,:);
    gab(:, :, i) = [Ri(1,1), Ri(1,2), Ri(1,3), di(1); ...
                    Ri(2,1), Ri(2,2), Ri(2,3), di(2); ... 
                    Ri(3,1), Ri(3,2), Ri(3,3), di(3); ...
                    0, 0, 0, 1];
end

% now we need to find g_opt,k
% extrapolate from initial position of marker
a0 = opti_change(4, 1);
b0 = opti_change(1, 1);
c0 = opti_change(2, 1);
d0 = opti_change(3, 1);
Ropt0 = [[a0*a0+b0*b0-c0*c0-d0*d0, 2*(b0*c0-a0*d0), 2*(b0*d0 + a0*c0)];...
         [2*(b0*c0+a0*d0), a0*a0-b0*b0+c0*c0-d0*d0, 2*(c0*d0-a0*b0)];...
         [2*(b0*d0-a0*c0), 2*(c0*d0+a0*b0), a0*a0-b0*b0-c0*c0+d0*d0]];
gopt0 = [Ropt0(1,1), Ropt0(1,2), Ropt0(1,3), rb_xyz(1,1); ...
         Ropt0(2,1), Ropt0(2,2), Ropt0(2,3), rb_xyz(2,1); ... 
         Ropt0(3,1), Ropt0(3,2), Ropt0(3,3), rb_xyz(3,1); ...
         0, 0, 0, 1];
     
% now compute estimate at frame k
gest = zeros(4,4,ndifferences + 1);
gest(:,:,1) = gopt0;
for i = 1:ndifferences
    gest(:,:,i+1) = gab(:,:,i)*gest(:,:,i);
end

disp = squeeze(gest(1:3, 4, :));


%% Plots
figure
subplot(2,2,1)
plot(opti_change(1,:), '--', 'LineWidth', 1.4)
hold on;
plot(q(1,:))
ylabel('q_0 (column W)')
xlabel('Frame Index')

subplot(2,2,2)
plot(opti_change(2,:), '--', 'LineWidth', 1.4)
hold on;
plot(q(2,:))
ylabel('q_1 (column X)')
xlabel('Frame Index')

subplot(2,2,3)
plot(opti_change(3,:), '--', 'LineWidth', 1.4)
hold on;
plot(q(3,:))
ylabel('q_2 (column Y)')
xlabel('Frame Index')

subplot(2,2,4)
plot(opti_change(4,:), '--', 'LineWidth', 1.4)
hold on;
plot(q(4,:))
ylabel('q_3 (column Z)')
xlabel('Frame Index')

legend('OptiTrack', 'Estimated')

sgtitle('Change in orientation between adjacent pairs of frames') 

figure
subplot(2,2,1)
plot(rb_q(1,:), '--', 'LineWidth', 1.4)
hold on;
plot(q2(1,:))
ylabel('q_0 (column W)')
xlabel('Frame Index')

subplot(2,2,2)
plot(rb_q(2,:), '--', 'LineWidth', 1.4)
hold on;
plot(q2(2,:))
ylabel('q_1 (column X)')
xlabel('Frame Index')

subplot(2,2,3)
plot(rb_q(3,:), '--', 'LineWidth', 1.4)
hold on;
plot(q2(3,:))
ylabel('q_2 (column Y)')
xlabel('Frame Index')

subplot(2,2,4)
plot(rb_q(4,:), '--', 'LineWidth', 1.4)
hold on;
plot(q2(4,:))
ylabel('q_3 (column Z)')
xlabel('Frame Index')

legend('OptiTrack', 'Estimated')

sgtitle('Change in orientation relative to first frame') 

% 2D
figure
rb_q(:,1200) = NaN; % cleanup
plot(time, rb_q,'--','LineWidth',1.4)
hold on;
plot(time, q2)
legend('q_w','q_x','q_y','q_z', 'q_w (Est.)','q_x (Est.)','q_y (Est.)','q_z (Est.)')
xlabel('Time(s)')
ylabel('Quaternion components')
title("Estimating Quaternion Components");


% Uncomment once disp = (x, y, z) is calculated

rb_xyz(:,1200) = NaN;
figure
plot(time,rb_xyz,'--','LineWidth',1.4)
hold on;
plot(time, disp)
xlabel('Time(s)')
ylabel('World Position')
legend('x','y','z','x (Est.)', 'y (Est.)', 'z (Est.)')

% 3D Path
figure
plot3(rb_xyz(1,:),rb_xyz(2,:),rb_xyz(3,:),'--','LineWidth',1.4)
hold on
plot3(disp(1,:),disp(2,:),disp(3,:))
axis equal
grid on
legend('Opti-Track', 'Calculated','Location','s')
title("Estimating the 3D Path");
