% Read OptiTrack File
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

%% Quaternion-Based Least Squares Approach
% Compute the centroids of the points
centroid = zeros(size(marker{1}));
for i=1:nmarkers
   centroid  = centroid + marker{i};
end
centroid = centroid./nmarkers;
Q_cent = centroid(2:end,:);
P_cent = centroid(1:end-1,:);

% Construct Q and P matrices
P = cell(nmarkers,1); % P'
Q = cell(nmarkers,1); % Q'
for i=1:nmarkers
   P{i} = marker{i}(1:end-1,:) - P_cent;
   Q{i} = marker{i}(2:end,:) - Q_cent;
end

Pmat = cell(nmarkers,1);
Qmat = cell(nmarkers,1);
npoints = size(P{1},1);
for i=1:nmarkers
    for j = 1:npoints
        Pix = P{i}(j,1) ;
        Piy = P{i}(j,2);
        Piz = P{i}(j,3);
        Pmat{i}(:,:,j) = [0   -Pix -Piy -Piz ;...
                          Pix  0    Piz -Piy;...
                          Piy -Piz  0    Pix ;...
                          Piz  Piy -Pix  0   ];
        Qix = Q{i}(j,1) ;
        Qiy = Q{i}(j,2);
        Qiz = Q{i}(j,3);
        Qmat{i}(:,:,j) = [0   -Qix -Qiy -Qiz;...
                          Qix  0   -Qiz  Qiy;...
                          Qiy  Qiz  0   -Qix;...
                          Qiz -Qiy  Qix  0];
    end
end

N = zeros(size(Pmat{1}));
for j = 1:npoints
    for i=1:nmarkers
        N(:,:,j) = N(:,:,j) + Pmat{i}(:,:,j) * transpose(Qmat{i}(:,:,j));
    end
end

% Find eigenvector of N associated with largest eigenvalue of N
q = zeros(4,npoints);
for j = 1:npoints
    [V, D] = eig(N(:,:,j));
    [d,ind] = sort(diag(D),'descend'); % Sort eigenvalues
    Ds = D(ind,ind); % Reorder eigenvalue matrix
    Vs = V(:,ind); % Reorder eigenvectors
    q(:,j) = Vs(:,1);
    if q(1,j) < 0 % Fix sign
        q(:,j) = -q(:,j);
    end
    q(:,j) = quatnormalize(q(:,j)'); % Normalize - unit quaternion
end
rb_q = [ rb_q(4,:); rb_q(1,:); rb_q(2,:); rb_q(3,:)];

% Find displacements between each adjacent pair of frames i
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


%% Plots
start = 1201;
stop = 2400;
figure
subplot(2,2,1)
plot(opti_change(1,start:stop), '--', 'LineWidth', 1.4)
hold on;
plot(q(1,start:stop))
ylabel('q_0 (column W)')
xlabel('Frame Index')
%legend('OptiTrack', 'Estimated')

subplot(2,2,2)
plot(opti_change(2,start:stop), '--', 'LineWidth', 1.4)
hold on;
plot(q(2,start:stop))
ylabel('q_1 (column X)')
xlabel('Frame Index')
%legend('OptiTrack', 'Estimated')

subplot(2,2,3)
plot(opti_change(3,start:stop), '--', 'LineWidth', 1.4)
hold on;
plot(q(3,start:stop))
ylabel('q_2 (column Y)')
xlabel('Frame Index')
%legend('OptiTrack', 'Estimated')

subplot(2,2,4)
plot(opti_change(4,start:stop), '--', 'LineWidth', 1.4)
hold on;
plot(q(4,start:stop))
ylabel('q_3 (column Z)')
xlabel('Frame Index')
%legend('OptiTrack', 'Estimated')

sgtitle('Change in orientation between adjacent pairs of frames') 

diff = max(sum(q(:, 1:1199) - opti_change(:, 1:1199)));

%%
start = 1201;
stop = 2400;
figure
subplot(2,2,1)
plot(rb_q(1,start:stop), '--', 'LineWidth', 1.4)
hold on;
plot(q2(1,start:stop))
ylabel('q_0 (column W)')
xlabel('Frame Index')
legend('OptiTrack', 'Estimated')

subplot(2,2,2)
plot(rb_q(2,start:stop), '--', 'LineWidth', 1.4)
hold on;
plot(q2(2,start:stop))
ylabel('q_1 (column X)')
xlabel('Frame Index')
legend('OptiTrack', 'Estimated')

subplot(2,2,3)
plot(rb_q(3,start:stop), '--', 'LineWidth', 1.4)
hold on;
plot(q2(3,start:stop))
ylabel('q_2 (column Y)')
xlabel('Frame Index')
legend('OptiTrack', 'Estimated')

subplot(2,2,4)
plot(rb_q(4,start:stop), '--', 'LineWidth', 1.4)
hold on;
plot(q2(4,start:stop))
ylabel('q_3 (column Z)')
xlabel('Frame Index')
legend('OptiTrack', 'Estimated')

sgtitle('Change in orientation relative to first frame') 

diff_first = max(sum(rb_q(:, 1:1199) - q2(:, 1:1199)));
%%

% 2D
figure
plot(time(start:stop), rb_q(:, start:stop))
hold on;
plot(time(start:stop), q2(:, start:stop), '--', 'LineWidth', 1.4)
legend('q_w','q_x','q_y','q_z', 'q_w (Est.)','q_x (Est.)','q_y (Est.)','q_z (Est.)')
xlabel('Time(s)')
ylabel('Quaternion components')

title('Estimated vs actual quaternion components') 

diff_quat = max(sum(q2(:, 1:1199) - rb_q(:, 1:1199)));

%% x,y,z
% Uncomment once disp = (x, y, z) is calculated 
%calculating xyz
d_12 = zeros(2401, 3);
for i = 1:1199
    phi = 2*acos(opti_change(4, i));
    w = [opti_change(1, i)/sin(phi/2), opti_change(2, i)/sin(phi/2), opti_change(3, i)/sin(phi/2)];
    x = w(1);
    y = w(2);
    z = w(3);
    a = opti_change(4, i);
    b = opti_change(1, i);
    c = opti_change(2, i);
    d = opti_change(3, i);
    R = [[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d + a*c)];...
        [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)];...
        [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]];
    R_12 = [[cos(phi) + x^2*(1-cos(phi)), x*y*(1-cos(phi)- z*sin(phi)), x*z*(1-cos(phi)) + y*sin(phi)];...
        [y*x*(1-cos(phi)) + z*sin(phi), cos(phi) + y^2*(1-cos(phi)), y*z*(1-cos(phi)) - x*sin(phi)];...
        [z*x*(1-cos(phi)) - y*sin(phi), z*y*(1-cos(phi)) + x*sin(phi), cos(phi) + z^2*(1-cos(phi))]];
    d_12(i, :) = Q_cent(i, :)' - R*P_cent(i, :)';
end

pos_first = zeros(2401, 3);
pos_first(1, :) = reshape(rb_xyz(:, 1),[1, 3]);
for i = 1:1199
    pos_first(i+1, :) = pos_first(i, :) + (d_12(i+1, :) - d_12(i, :));
end

% rb_trans = zeros(3, 2401);
% for i = 2:1200
%     rb_trans(:, i-1) = rb_xyz(:, i) - rb_xyz(:, i-1);
% end

disp = d_12;
figure
 plot(time,rb_xyz)
 hold on;
 plot(time, disp)
 xlabel('Time(s)')
 ylabel('World Position')
 legend('x','y','z', 'x_{est}', 'y_{est}', 'z_{est}')
% 
 % 3D Path
 figure
 plot3(rb_xyz(1,:),rb_xyz(2,:),rb_xyz(3,:))
 hold on
 plot3(pos_first(1:1199,1),pos_first(1:1199,2),pos_first(1:1199,3))
 axis equal
 grid on


