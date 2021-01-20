close all
clear
clc


##############################################################################################
################################### LOAD SUPPORT LIBRARIES ###################################
##############################################################################################

addpath './g2o_wrapper'
source './utils.m'
source './poses_linear_system.m'
source './bearings_linear_system.m'
source './bearing_only_slam_least_square.m'


##############################################################################################
##################################### LOAD GROUND TRUTH ######################################
##############################################################################################

# load ground truth dataset

[landmarks, poses, transitions, observations] = loadG2o('slam2D_bearing_only_ground_truth.g2o');


# there are some problems with load2Go compatibility with new octave version,
# first positions of each struct are void, thus cut and take the rest

landmarks = landmarks(2:end);
poses = poses(2:end);
transitions = transitions(2:end);
observations = observations(2:end);


# define some global variables

global num_poses = size(poses,2);
global num_landmarks = size(landmarks,2);
global pose_dim = 3;
global landmark_dim = 2;
global system_size = num_poses*pose_dim+num_landmarks*landmark_dim;
global num_observations = size(observations,2);


# the state is made by X={XR1,..XRN,XL1,..XLM}
# XR is in SE(2), take each (x y theta) vector and map it into [R|t]
# then map the indices with ids

global XR_true = zeros(3,3,num_poses);
global XL_true = zeros(2,num_landmarks);

# id_to_index maps the pose or landmark id to it's respective index

global xr_id_to_index = zeros(1,num_poses);
global xl_id_to_index = zeros(1,num_landmarks);

# collect poses ground truth

for i=1:num_poses
	global XR_true;
	global xr_id_to_index;
	XR_true(:,:,i)=v2t([poses(i).x poses(i).y poses(i).theta]);
	xr_id_to_index(:,i)=poses(i).id;
endfor;


# XL is in R2, just copy and map the indices with the ids

# collect landmarks ground truth

for i=1:num_landmarks
	global XL_true;
	global xl_id_to_index;
	XL_true(:,i)=[landmarks(i).x_pose landmarks(i).y_pose];
	xl_id_to_index(:,i)=landmarks(i).id;
endfor;


##############################################################################################
##################################### LOAD INITIAL GUESS #####################################
##############################################################################################


global XR_guess = zeros(3,3,num_poses);
global XL_guess = zeros(2,num_landmarks);

# load initial guess dataset, take initial pose guess and initialize the landmarks by
# parsing the observations

[landmarks, poses, transitions, observations] = loadG2o('slam2D_bearing_only_initial_guess.g2o');


# same compatibility problem as before, cut first void element

landmarks = landmarks(2:end);
poses = poses(2:end);
transitions = transitions(2:end);
observations = observations(2:end);


# collect initial guess

for i=1:num_poses
	global XR_guess;
	global xr_id_to_index;
	XR_guess(:,:,i)=v2t([poses(i).x poses(i).y poses(i).theta]);
	xr_id_to_index(:,i)=poses(i).id;
endfor;


##############################################################################################
################################### PARSE LANDMARKS EDGES ####################################
##############################################################################################

# count total observations for parse and initializing flattened observations vector

global num_measurements = 0;
for i=1:num_observations
	num_observed = size(observations(i).observation,2);
	for j=1:num_observed
		num_measurements++; # total measurements
	endfor;
endfor;


# flattened observations is made by three rows (pose_index, land_index, bearing) and num_measurements cols

global ZL = zeros(1,num_measurements);
global associations_ZL = zeros(2,num_measurements);

num_measurement = 1;
for i=1:num_observations
	observed = observations(i);
	pose_id = observed.pose_id;
	pose_index = index_from_id(xr_id_to_index,pose_id);			# get index corresponding to id
	num_observed = size(observed.observation,2);
	for j=1:num_observed
		measurement = observed.observation(j);
		landmark_id = measurement.id;
		landmark_index = index_from_id(xl_id_to_index,landmark_id);	# get index corresponding to id
		ZL(1,num_measurement)= measurement.bearing;
		associations_ZL(:,num_measurement) = [pose_index landmark_index]';
		num_measurement++;
	endfor;
endfor;


##############################################################################################
##################################### PARSE POSES EDGES ######################################
##############################################################################################

# count total transitionss

global num_transitions = size(transitions,2);

# flattened transitions is made by 3x3 homogeneous matrices and associations are made by (pose_index_i, pose_index_j)

global ZR = zeros(3,3,num_transitions);
global assosiations_ZR = zeros(2,num_transitions);

for i=1:num_transitions
	pose_i_id = transitions(i).id_from;
	pose_j_id = transitions(i).id_to;
	pose_i_index = index_from_id(xr_id_to_index,pose_i_id);			# get index corresponding to id
	pose_j_index = index_from_id(xr_id_to_index,pose_j_id);
	associations_ZR(:,i)=[pose_i_index pose_j_index]';
	ZR(:,:,i)=v2t(transitions(i).v);
endfor;


##############################################################################################
################################ LANDMARKS LINEAR TRIANGULATION ##############################
##############################################################################################

# for each landmark store the relative robot pose and bearings
# indices : ZL indices sorted by landmarks
# out : sorted landmarks, repetitions are the number of observations per landmarks
# landmarks_ids_list : list of landmarks ids
# observations_per_landmarks : store how many observations per each landmark

[out,indices]=sort(associations_ZL(2,:));
landmarks_ids_list = unique(out);
observations_per_landmarks = arrayfun(@(x)sum(out==x), landmarks_ids_list);

# for each landmark get the relative bearing measurements and consistent poses and find the point of intersection
# of the lines passing througt the two most orthogonal poses w.r.t the relative landmark

poses_relative_distances_threshold = 5; # distance treshold which cutt off potential outliers (need to be tuned)

for i=1:size(landmarks_ids_list,2)
	landmark_id =landmarks_ids_list(:,i);
	num_measurements = observations_per_landmarks(:,i);
	offset = sum(observations_per_landmarks(:,1:i))-num_measurements;
	distances = zeros(num_measurements);
	positions = zeros(3, num_measurements);
	bearings = zeros(1,num_measurements);
	for j=1:num_measurements
		index = indices(:,offset+j);
		positions(:,j) = t2v(XR_guess(:,:,associations_ZL(1,index)));
		bearings(:,j) = ZL(1,index);
	endfor;
	# get matrix of pose relative distances to cut off outliers given a threshold
	# we can get matrix since distance is independent is the order of factors
	if num_measurements > 1
		dist_matrix = sqrt((positions(1,:)-positions(1,:)').^2+(positions(2,:)-positions(2,:)').^2);
		outliers_index = find(sum(dist_matrix>poses_relative_distances_threshold,1)==num_measurements-1);
		positions(:,outliers_index)=[];					# remove outliers
		bearings(:,outliers_index)=[];
		num_measurements-=size(outliers_index,2);
	endif;
	positions(3,:) = bearings(:,:)+positions(3,:);				# compute total angle of bearing w.r.t x axis
	positions(3,:) = atan2(sin(positions(3,:)),cos(positions(3,:)));	# normalize the sum in -pi and pi
	if num_measurements == 1
		XL_guess(:,i)= positions(1:2,1) + [cos(positions(3,1)) sin(positions(3,1))]';
	else
		diff_matrix = positions(3,:)-positions(3,:)';			# compute difference of each total angle w.r.t the others
		diff_matrix = atan2(sin(diff_matrix),cos(diff_matrix));		# normalize the difference in -pi and pi
		# compute the abs of the sin of the difference, this means the more the angles are orthogonal the higher the value
		sin_matrix = abs(sin(diff_matrix));
		sin_matrix = triu(sin_matrix,1);				# get the upper triangular matrix since we have done abs value
		[max_rows, max_cols] = find(sin_matrix==max(max(sin_matrix)));	# find the most orthogonal poses
		max_rows = max_rows(1);						# avoid cases in which there are more maximums
		max_cols = max_cols(1);
		p_a = positions(1:2,max_rows);
		a_a = positions(3,max_rows);
		p_aa = p_a + [cos(a_a) sin(a_a)]';				# compute second point
		p_b = positions(1:2,max_cols);
		a_b = positions(3,max_cols);
		p_bb = p_b + [cos(a_b) sin(a_b)]';				# compute second point
		XL_guess(:,i) = get2dLinesIntersection(p_a,p_aa,p_b,p_bb);	# compute intersection point
	endif;
endfor;


##############################################################################################
################################## RUN TOTAL LEAST SQUARE ####################################
##############################################################################################

# define optimization variables

global kernel_threshold = 10;
global damping = 1;
global num_iterations = 20;

disp('BearingOnlySLAMLeastSquare running...');

[XR, XL, chi_stats_l, inliers_stats_l, chi_stats_r, inliers_stats_r, H, b, dx] = bearingOnlySLAM(XR_guess,XL_guess,ZL,associations_ZL,ZR,associations_ZR,num_iterations,kernel_threshold,damping);

# plot results

plotState(XR_true, XR_guess, XR, XL_true, XL_guess, XL, chi_stats_l, inliers_stats_l, chi_stats_r, inliers_stats_r, H);
