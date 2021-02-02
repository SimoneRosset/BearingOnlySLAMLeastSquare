close all
clear
clc


##############################################################################################
################################### LOAD SUPPORT LIBRARIES ###################################
##############################################################################################

source './utils.m'
source './poses_linear_system.m'
source './bearings_linear_system.m'
source './bearing_only_slam_least_square.m'


##############################################################################################
###################################### TRANSITION MODEL ######################################
##############################################################################################

function xytheta_prime = transition_model(xytheta, uxutheta)

        xytheta_prime = xytheta;
        x = xytheta(1);
        y = xytheta(2);
	theta = xytheta(3);
        ux = uxutheta(1);
        utheta = uxutheta(2);
        c = cos(theta);
        s = sin(theta);

        %       (x' ) = (    x + ux*cos(th)     )
        %       (y' ) = (    y + ux*sin(th)     )
        %       (th') = (    th + uth           )

        xytheta_prime(1) = x + ux*c;
        xytheta_prime(2) = y + ux*s;
        xytheta_prime(3) = theta + utheta;
	xytheta_prime(3) = atan2(sin(xytheta_prime(3)),cos(xytheta_prime(3)));
endfunction;

##############################################################################################
#################################### CREATE GROUND TRUTH #####################################
##############################################################################################

# define some global variables

global num_poses = 301;
global num_landmarks = 141;
global pose_dim = 3;
global landmark_dim = 2;
global system_size = num_poses*pose_dim+num_landmarks*landmark_dim;
global num_observations = 300;
global num_transitions = 300;

# the state is made by X={XR1,..XRN,XL1,..XLM}
# XR is in SE(2), take each (x y theta) vector and map it into [R|t]
# then map the indices with ids

global XR_true = zeros(3,3,num_poses);
global XL_true = zeros(2,num_landmarks);
global poses = zeros(3,num_poses);
global landmarks = zeros(2,num_landmarks);

# id_to_index maps the pose or landmark id to it's respective index

global xr_id_to_index = zeros(1,num_poses);
global xl_id_to_index = zeros(1,num_landmarks);


##############################################################################################
###################################### POSE GROUND TRUTH #####################################
##############################################################################################

ux = [-1,0,1];
utheta = [0,pi/2,-pi/2];
XR_true(:,:,1)=eye(3);

for i=2:num_poses
	idx = randperm(length(ux),1);
	dux = ux(idx);
	idx = randperm(length(utheta),1);
	dtheta = utheta(idx);
	uxutheta =[dux,dtheta];
	xytheta = poses(:,i-1);
	xytheta_prime = transition_model(xytheta,uxutheta);
	poses(:,i)=xytheta_prime;
	XR_true(:,:,i)=v2t(xytheta_prime);
	xr_id_to_index(:,i)=i;
endfor;


##############################################################################################
###################################### POSE-POSE EDGES #######################################
##############################################################################################

global associations_ZR = [[i=1:num_transitions];[i=2:num_transitions+1]];

global ZR = zeros(3,3,num_transitions);

for i=1:num_transitions
	ZR(:,:,i)=XR_true(:,:,i)\XR_true(:,:,i+1);
endfor;


##############################################################################################
################################# LANDMARKS GROUND TRUTH #####################################
##############################################################################################

minx = min(poses(1,:));
miny = min(poses(2,:));
maxx = max(poses(1,:));
maxy = max(poses(2,:));
width = maxx-minx;
height = maxy-miny;

for i=1:num_landmarks
	landmarks(:,i)=[minx+rand(1)*width; miny+rand(1)*height];
	xl_id_to_index(:,i)=i;
endfor;

XL_true = landmarks;


##############################################################################################
#################################### POSE-LANDMARKS EDGES ####################################
##############################################################################################

dist_matrix = sqrt((poses(1,:)'-landmarks(1,:)).^2+(poses(2,:)'-landmarks(2,:)).^2);

bearing_range = 1;

[pose_i,land_j]=find(dist_matrix<bearing_range);

global associations_ZL = [pose_i'; land_j'];

[s,i] = sort(associations_ZL(1,:));

associations_ZL = associations_ZL(:,i);

global num_measurements = size(associations_ZL,2);

global ZL = zeros(1,num_measurements);

# XL is in R2, just copy and map the indices with the ids

# collect landmarks ground truth

for i=1:num_measurements
	pose_id = associations_ZL(1,i);
	land_id = associations_ZL(2,i);
	X = XR_true(:,:,pose_id);
	l = XL_true(:,land_id);
	xl = X\[l;1];
	ZL(:,i)=atan2(xl(2),xl(1));
endfor;

##############################################################################################
################################## RUN TOTAL LEAST SQUARE ####################################
##############################################################################################

# define optimization variables

global kernel_threshold = 10;
global damping = 1;
global num_iterations = 10;

disp('BearingOnlySLAMLeastSquare running...');

XR_guess = XR_true;
XL_guess = XL_true;

[XR, XL, chi_stats_l, inliers_stats_l, chi_stats_r, inliers_stats_r, H, b, dx] = bearingOnlySLAM(XR_guess,XL_guess,ZL,associations_ZL,ZR,associations_ZR,num_iterations,kernel_threshold,damping);

# plot results


function plotState(XR_true, XR_guess, XR, XL_true, XL_guess, XL, chi_stats_l, inliers_stats_l, chi_stats_r, inliers_stats_r, H)
        # Plot State
        figure(1);
        hold on;
        grid;

        subplot(2,2,1);
        title("Landmark Initial Guess");
        plot(XL_true(1,:),XL_true(2,:),'b*',"linewidth",2);
        hold on;
        plot(XL_guess(1,:),XL_guess(2,:),'ro',"linewidth",2);
        legend("Landmark True", "Guess");grid;


        subplot(2,2,2);
        title("Landmark After Optimization");
        plot(XL_true(1,:),XL_true(2,:),'b*',"linewidth",2);
        hold on;
        plot(XL(1,:),XL(2,:),'ro',"linewidth",2);
        legend("Landmark True", "Guess");grid;


        subplot(2,2,3);
        title("Poses Initial Guess");
        plot(squeeze(XR_true(1,3,:)),squeeze(XR_true(2,3,:)),'b*-',"linewidth",2);
        hold on;
        plot(squeeze(XR_guess(1,3,:)),squeeze(XR_guess(2,3,:)),'ro-',"linewidth",2);
        legend("Poses True", "Guess");grid;


        subplot(2,2,4);
        title("Poses After Optimization");
        plot(squeeze(XR_true(1,3,:)),squeeze(XR_true(2,3,:)),'b*-',"linewidth",2);
        hold on;
        plot(squeeze(XR(1,3,:)),squeeze(XR(2,3,:)),'ro-',"linewidth",2);
        legend("Poses True", "Guess"); grid;



        figure(2);
        hold on;
        grid;
        title("chi evolution");

        subplot(3,2,1);
        plot(chi_stats_r, 'r-', "linewidth", 2);
        legend("Chi Poses"); grid; xlabel("iterations");
        subplot(3,2,2);
        plot(inliers_stats_r, 'b-', "linewidth", 2);
        legend("#poses inliers"); grid; xlabel("iterations");

        subplot(3,2,3);
        plot(chi_stats_l, 'r-', "linewidth", 2);
        legend("Chi Landmark"); grid; xlabel("iterations");
        subplot(3,2,4);
        plot(inliers_stats_l, 'b-', "linewidth", 2);
        legend("#landmarks inliers"); grid; xlabel("iterations");

        subplot(3,2,5);
        plot(chi_stats_l+chi_stats_r, 'r-', "linewidth", 2);
        legend("Chi total"); grid; xlabel("iterations");
        subplot(3,2,6);
        plot(inliers_stats_l+inliers_stats_r, 'b-', "linewidth", 2);
        legend("#total inliers"); grid; xlabel("iterations");

        figure(3);
        title("H matrix");
        H_ =  H./H;                      # NaN and 1 element
        H_(isnan(H_))=0;                 # Nan to Zero
        H_ = abs(ones(size(H_)) - H_);   # switch zero and one
        H_ = flipud(H_);                 # switch rows
        colormap(gray(64));
        hold on;
        image([0.5, size(H_,2)-0.5], [0.5, size(H_,1)-0.5], H_*64);
        hold off;
endfunction;

plotState(XR_true, XR_guess, XR, XL_true, XL_guess, XL, chi_stats_l, inliers_stats_l, chi_stats_r, inliers_stats_r, H);
