1;

##############################################################################################
#################################### LEAST SQUARE POSES ######################################
##############################################################################################


# error and jacobian of a measured pose, all poses are in world frame
# input:
#   Xi: the observing robot pose (3x3 homogeneous matrix)
#   Xj: the observed robot pose (3x3 homogeneous matrix)
#   Z:   the relative transform measured between Xr1 and Xr2
#   e: 6x1 is the difference between prediction, and measurement, vectorized
#   Ji : 6x3 derivative w.r.t a the error and a perturbation of the
#       first pose
#   Jj : 6x3 derivative w.r.t a the error and a perturbation of the
#       second pose

function [e,Ji,Jj]=poseErrorAndJacobian(Xi,Xj,Z)
        Ri=Xi(1:2,1:2);
        Rj=Xj(1:2,1:2);
        ti=Xi(1:2,3);
        tj=Xj(1:2,3);
        R_der_0=[0 -1;1 0];

        Z_hat=eye(3);
        Z_hat(1:2,1:2)=Ri'*Rj;
        Z_hat(1:2,3)=Ri'*(tj-ti);
        e=flattenIsometryByColumns(Z_hat-Z);

        Ji=zeros(6,3);
        Jj=zeros(6,3);

        Jj(1:4,3)=reshape(Ri'*R_der_0*Rj, 4, 1);
        Jj(5:6,1:2)=Ri';

        Jj(5:6,3)=Ri'*R_der_0*tj;
        Ji=-Jj;
endfunction;

#linearizes the robot-robot measurements
# inputs:
#   XR: the initial robot poses (3x3xnum_poses: array of homogeneous matrices)
#   XL: the initial landmark estimates (2xnum_landmarks matrix of landmarks)
#   ZR: the robot_robot measuremenrs (3x3xnum_measurements: array of homogeneous matrices)
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   kernel_threshod: robust kernel threshold
# outputs:
#   H: the H matrix, filled
#   b: the b vector, filled
#   chi_tot: the total chi2 of the current round
#   num_inliers: number of measurements whose error is below kernel_threshold

function [H,b, chi_tot, num_inliers]=buildLinearSystemPoses(XR, XL, ZR, associations_ZR, kernel_threshold)
        global pose_dim;
        global landmark_dim;
        num_poses = size(XR, 3);
        num_landmarks = size(XL,2);
        system_size=pose_dim*num_poses+landmark_dim*num_landmarks;
        H=zeros(system_size, system_size);
        b=zeros(system_size,1);
        chi_tot=0;
        num_inliers=0;
        num_measurements = size(ZR,3);

        for (i=1:num_measurements)
                Omega=eye(6);
                Omega(1:4,1:4)*=1e-5;
                pose_i_index=associations_ZR(1,i);
                pose_j_index=associations_ZR(2,i);
                z=ZR(:,:,i);
                Xi=XR(:,:,pose_i_index);
                Xj=XR(:,:,pose_j_index);
                [e,Ji,Jj] = poseErrorAndJacobian(Xi, Xj, z);
                chi=e'*Omega*e;
                if (chi>kernel_threshold)
                        Omega*=sqrt(kernel_threshold/chi);
                        chi=kernel_threshold;
                else
                        num_inliers ++;
                endif;
                chi_tot+=chi;

                pose_i_matrix_index=poseMatrixIndex(pose_i_index, num_poses, num_landmarks);
                pose_j_matrix_index=poseMatrixIndex(pose_j_index, num_poses, num_landmarks);

                H(pose_i_matrix_index:pose_i_matrix_index+pose_dim-1,
                pose_i_matrix_index:pose_i_matrix_index+pose_dim-1)+=Ji'*Omega*Ji;

                H(pose_i_matrix_index:pose_i_matrix_index+pose_dim-1,
                pose_j_matrix_index:pose_j_matrix_index+pose_dim-1)+=Ji'*Omega*Jj;

                H(pose_j_matrix_index:pose_j_matrix_index+pose_dim-1,
                pose_i_matrix_index:pose_i_matrix_index+pose_dim-1)+=Jj'*Omega*Ji;

                H(pose_j_matrix_index:pose_j_matrix_index+pose_dim-1,
                pose_j_matrix_index:pose_j_matrix_index+pose_dim-1)+=Jj'*Omega*Jj;

                b(pose_i_matrix_index:pose_i_matrix_index+pose_dim-1)+=Ji'*Omega*e;
                b(pose_j_matrix_index:pose_j_matrix_index+pose_dim-1)+=Jj'*Omega*e;
        endfor;
endfunction;
