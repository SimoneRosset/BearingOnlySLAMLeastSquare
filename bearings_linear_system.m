1;

##############################################################################################
################################### LEAST SQUARE BEARING #####################################
##############################################################################################


# error and jacobian of a measured landmark
# input:
#   Xr: the robot pose in world frame (3x3 homogeneous matrix)
#   Xl: the landmark pose (2x1 vector, 2d pose in world frame)
#   z:  measured bearing of landmark
# output:
#   e: 1x1 is the difference between prediction and measurement
#   Jr: 1x3 derivative w.r.t a the error and a perturbation on the
#       pose
#   Jl: 1x2 derivative w.r.t a the error and a perturbation on the
#       landmark

function [e,Jr, Jl]=bearingErrorAndJacobian(Xr,Xl,z)
        t=Xr(1:2,3);
        R=Xr(1:2,1:2);
        R_der_0=[0 -1; 1 0];
        p_hat=R'*(Xl-t);
        z_hat=atan2(p_hat(2),p_hat(1));
        e=z_hat-z;
        e=atan2(sin(e), cos(e)); # normalize error
        J_r=zeros(2,3);
        J_l=zeros(2,2);
        Jr=zeros(1,3);
        Jl=zeros(1,2);
        J_r(1:2,1:2)=-R';               # derivative w.r.t Delta t (position increment)
        J_r(1:2,3)=R'*R_der_0'*Xl;      # derivative w.r.t Delta theta
        J_l(:,:)=R';                    # derivative w.r.t Delta l (landmark increment)
        J_a=J_atan2(p_hat);             # atan derivative
        Jr=J_a*J_r;                     # chain rule
        Jl=J_a*J_l;                     # chain rule
endfunction;

#linearizes the robot-landmark measurements
#   XR: the initial robot poses (3x3xnum_poses: array of homogeneous matrices)
#   XL: the initial landmark estimates (2xnum_landmarks matrix of landmarks)
#   Z:  the measurements (1xnum_measurements)
#   kernel_threshod: robust kernel threshold
# output:
#   XR: the robot poses after optimization
#   XL: the landmarks after optimization
#   chi_tot: array 1:num_iterations, containing evolution of chi2
#   num_inliers: array 1:num_iterations, containing evolution of inliers
#   H: H matrix (system_sizexsystem_size) of linear system
#   b: b vector (system_sizex1)

function [H, b, chi_tot, num_inliers]=buildLinearSystemLandmarks(XR,XL,ZL,associations_ZL,kernel_threshold)
        global system_size;
        global pose_dim;
        global landmark_dim;
        num_poses = size(XR,3);
        num_landmarks = size(XL,2);
        num_measurements = size(ZL,2);
        H=zeros(system_size,system_size);
        b=zeros(system_size,1);
        chi_tot=0;
        chi=0;
        num_inliers=0;
        for i=1:num_measurements
                pose_index = associations_ZL(1,i);
                landmark_index = associations_ZL(2,i);
                Xr = XR(:,:,pose_index);
                Xl = XL(:,landmark_index);
                z = ZL(1,i);
                [e,Jr, Jl]=bearingErrorAndJacobian(Xr,Xl,z);
                chi=e'*e;
                if (chi>kernel_threshold)
                        e*=sqrt(kernel_threshold/chi);
                        chi=kernel_threshold;
                else
                        num_inliers++;
                endif;
                chi_tot+=chi;
                pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
                landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);

                H(pose_matrix_index:pose_matrix_index+pose_dim-1,
                pose_matrix_index:pose_matrix_index+pose_dim-1)+=Jr'*Jr;

                H(pose_matrix_index:pose_matrix_index+pose_dim-1,
                landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Jr'*Jl;

                H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
                landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Jl'*Jl;

                H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
                pose_matrix_index:pose_matrix_index+pose_dim-1)+=Jl'*Jr;

                b(pose_matrix_index:pose_matrix_index+pose_dim-1)+=Jr'*e;
                b(landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Jl'*e;
        endfor;
endfunction;
