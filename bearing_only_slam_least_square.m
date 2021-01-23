1;


##############################################################################################
############################ LEAST SQUARE BEARING ONLY SLAM  #################################
##############################################################################################



# implementation of the boxplus
# applies a perturbation to a set of landmarks and robot poses
# input:
#   XR: the robot poses (4x4xnum_poses: array of homogeneous matrices)
#   XL: the landmark pose (3xnum_landmarks matrix of landmarks)
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   dx: the perturbation vector of appropriate dimensions
#       the poses come first, then the landmarks
# output:
#   XR: the robot poses obtained by applying the perturbation
#   XL: the landmarks obtained by applying the perturbation

function [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx)
        global pose_dim;
        global landmark_dim;
        for(pose_index=1:num_poses)
                pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
                dxr=dx(pose_matrix_index:pose_matrix_index+pose_dim-1,:);
                XR(:,:,pose_index)=v2t(dxr)*XR(:,:,pose_index);
        endfor;
        for(landmark_index=1:num_landmarks)
                landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);
                dxl=dx(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,:);
                XL(:,landmark_index)+=dxl;
        endfor;
endfunction;


# implementation of the optimization loop with robust kernel
# applies a perturbation to a set of landmarks and robot poses
# input:
#   XR: the initial robot poses (3x3xnum_poses: array of homogeneous matrices)
#   XL: the initial landmark estimates (2xnum_landmarks matrix of landmarks)
#   Z:  the measurements (1xnum_measurements)
#   num_iterations: the number of iterations of least squares
#   damping:      damping factor (in case system not spd)
#   kernel_threshod: robust kernel threshold
# output:
#   XR: the robot poses after optimization
#   XL: the landmarks after optimization
#   chi_stats: array 1:num_iterations, containing evolution of chi2 for landmarks, projections and poses
#   inliers_stats: array 1:num_iterations, containing evolution of inliers landmarks, projections and poses

function [XR, XL, chi_stats_l, inliers_stats_l, chi_stats_r, inliers_stats_r, H, b, dx] = bearingOnlySLAM(XR,XL,ZL,associations_ZL,ZR,associations_ZR,num_iterations,kernel_threshold,damping)
        global pose_dim;
        chi_stats_l = chi_stats_r = zeros(1,num_iterations);
        inliers_stats_l = inliers_stats_r = zeros(1,num_iterations);
        global num_poses;
        global num_landmarks;
        global system_size;
        for i=1:num_iterations
                [H_l, b_l, chi_tot_l, num_inliers_l]=buildLinearSystemLandmarks(XR,XL,ZL,associations_ZL,kernel_threshold);
                chi_stats_l(:,i)=chi_tot_l;
                inliers_stats_l(:,i)=num_inliers_l;

                [H_r, b_r, chi_tot_r, num_inliers_r]=buildLinearSystemPoses(XR, XL, ZR, associations_ZR, kernel_threshold);
                chi_stats_r(:,i)=chi_tot_r;
                inliers_stats_r(:,i)=num_inliers_r;

                b=b_r+b_l;
                H=H_r+H_l+eye(system_size,system_size)*damping;
                dx=zeros(system_size,1);
		H((num_poses-1)*pose_dim+1:num_poses*pose_dim,:)=[];
                H(:,(num_poses-1)*pose_dim+1:num_poses*pose_dim)=[];
                b((num_poses-1)*pose_dim+1:num_poses*pose_dim)=[];
                dx=-H\b;			                                                                    # solve the linear system
		dx = [dx(1:(num_poses-1)*pose_dim)' zeros(1,pose_dim) dx((num_poses-1)*pose_dim+1:end)']';          # fix last dimension
                [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx);                         # apply update
        endfor;
endfunction;
