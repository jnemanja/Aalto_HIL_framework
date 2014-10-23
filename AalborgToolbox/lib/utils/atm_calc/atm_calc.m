% satellite polygons  --  [[x y z]' [x y z]' ...]
translate = @(pol, v) pol+v'*ones(1, size(pol, 2));
rotate = @(pol, r) angle2dcm(deg2rad(r(1)), ...
                    deg2rad(r(2)), deg2rad(r(3)), 'XYZ')*pol;

small_plate = [0.0 -0.05 -0.05;
               0.0  0.05 -0.05;
               0.0  0.05  0.05;
               0.0 -0.05  0.05]';

large_plate = [0.1035 0.0 -0.05;
              -0.1035 0.0 -0.05;
              -0.1035 0.0  0.05;
               0.1035 0.0  0.05]';
            
antenna = [0.0    0.0075 0.0;
           0.0   -0.0075 0.0;
           -0.15 -0.0075 0.0;
           -0.15  0.0075 0.0]';

plate_front    = translate(rotate(small_plate, [  0   0   0]), [ 0.1035 0 0]);
plate_back     = translate(rotate(small_plate, [  0 180   0]), [-0.1035 0 0]);
plate_left     = translate(rotate(large_plate, [  0   0 180]), [0 -0.05 0]);
plate_right    = translate(rotate(large_plate, [  0   0   0]), [0  0.05 0]);
plate_top      = translate(rotate(large_plate, [ 90   0   0]), [0  0 -0.05]);
plate_bottom   = translate(rotate(large_plate, [-90   0   0]), [0  0  0.05]);
antenna_top_f    = translate(rotate(antenna, [  0 -45 0]),  [-0.1035  0     0.05]);
antenna_bottom_f = translate(rotate(antenna, [180  45 0]),  [-0.1035  0    -0.05]);
antenna_left_f   = translate(rotate(antenna, [-90  0 -45]), [-0.1035 -0.05  0]);
antenna_right_f  = translate(rotate(antenna, [ 90  0  45]), [-0.1035  0.05  0]);
antenna_top_b    = translate(rotate(antenna, [180 -45 0]),  [-0.1035  0     0.05]);
antenna_bottom_b = translate(rotate(antenna, [  0  45 0]),  [-0.1035  0    -0.05]);
antenna_left_b   = translate(rotate(antenna, [ 90  0 -45]), [-0.1035 -0.05  0]);
antenna_right_b  = translate(rotate(antenna, [-90  0  45]), [-0.1035  0.05  0]);


satellite = {plate_front, plate_back, plate_left, plate_right, ...
                plate_top, plate_bottom, antenna_top_f, antenna_bottom_f,...
                antenna_left_f, antenna_right_f, antenna_top_b,...
                antenna_bottom_b, antenna_left_b, antenna_right_b};
           
sat_rot = [0 0 40];

satellite = cellfun(@(pol) rotate(pol, sat_rot), satellite, 'UniformOutput', false);

mass_center = rotate([0.01 0 0]', sat_rot);

% plot figure
figure;

% plot satellite
subplot(2, 2, 1);
xlabel('x'); ylabel('y'); zlabel('z');  
axis([-0.3 0.2 -0.2 0.2 -0.2 0.2])
grid
hold on
view(60, 50)
cellfun(@(pol) draw_polygon(pol, 1), satellite);
plot3(mass_center(1), mass_center(2), mass_center(3), '*');
alpha(0.7);


%% back-face culling

satellite = cellfun(@(pol) polygon_backface_check(pol),...
    satellite, 'UniformOutput', false);
satellite = flatten(satellite);


%% polygon division

% % center of mass division
% satellite = cellfun(@(pol) polygon_plane_division(pol, [0 1 0]', mass_center), ...
%     satellite, 'UniformOutput', false);
% satellite = flatten(satellite);
% satellite = cellfun(@(pol) polygon_plane_division(pol, [0 0 1]', mass_center), ...
%     satellite, 'UniformOutput', false);
% satellite = flatten(satellite);

% % occulsion divison
% x = [1 0 0];
% for i = length(satellite)
%     pol = satellite{i};
%     lines = [pol; pol(:,2:end) pol(:,1)];
%     for line = lines(:, 1:end)
%         line_n = line(1:3) - line(4:6);
%         line_n = line_n/norm(line_n);
%         plane_n = cross(line_n, x);
%         satellite = cellfun(@(pol) polygon_plane_division(pol, plane_n, line(1:3)), ...
%             satellite, 'UniformOutput', false);
%         satellite = flatten(satellite);
%     end
% end

% plot divided satellite
subplot(2, 2, 3);
xlabel('x'); ylabel('y'); zlabel('z');  
axis([-0.3 0.2 -0.2 0.2 -0.2 0.2])
grid
hold on
view(90, 0)
cellfun(@(pol) draw_polygon(pol, 1), satellite);

%% mass to pressure vectors

satellite_c = cellfun(@(pol) polygon3_centroid(pol), satellite, 'UniformOutput', false);
rmp = cellfun(@(c) c'-mass_center, satellite_c, 'UniformOutput', false);

% plot centroids
cellfun(@(pol) plot3(pol(1), pol(2), pol(3), '+'), satellite_c);
alpha(0.7);


%% projection
 
% projection matrix
nv = [1 0 0];
norm_v = nv/norm(nv);
s = sum(norm_v.^2);
% proj_mat = [norm_v(2)^2+norm_v(3) -norm_v(1)*norm_v(2)   -norm_v(1)*norm_v(3);
%            -norm_v(1)*norm_v(2)    norm_v(1)^2+norm_v(3) -norm_v(2)*norm_v(3);
%            -norm_v(1)*norm_v(3)   -norm_v(2)*norm_v(3)    norm_v(1)^2+norm_v(2);]/s;
proj_mat = [0 0 0;
            0 1 0;
            0 0 1];
satellite_projected = cellfun(@(pol) proj_mat*pol, satellite, 'UniformOutput', false);

%% polygon areas

satellite_A = cellfun(@(pol) polygon2_area(pol(2:3,:)), satellite_projected,...
    'UniformOutput', false);
satellite_T = cellfun(@(r, A) cross(r, [-1*A*200 0 0]'),...
    rmp, satellite_A, 'UniformOutput', false);
 resultant_T = sum([satellite_T{:}],2);


% plot projected satellite
subplot(2, 2, 2);
xlabel('x'); ylabel('y'); zlabel('z');  
axis([-0.3 0.2 -0.2 0.2 -0.2 0.2])
grid
hold on
view(90, 0)
cellfun(@(pol, A) draw_polygon(pol, 5*A), satellite_projected, satellite_A);
subplot(2, 2, 3);
cellfun(@(T, C) line([C(1); C(1)+T(1)], [C(2); C(2)+T(2)], [C(3); C(3)+T(3)]),...
    satellite_T, satellite_c);
 
line([mass_center(1) mass_center(1)+resultant_T(1)],...
     [mass_center(2) mass_center(2)+resultant_T(2)],...
     [mass_center(3) mass_center(3)+resultant_T(3)],...
     'Color', 'red');


