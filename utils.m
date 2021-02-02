1;

##############################################################################################
######################################### UTILITIES ##########################################
##############################################################################################

# given an array of ids return the id query relative index
# inputs:
#   id_to_index : array of ids
#   id : query id
# output:
#   index : id relative index

function index = index_from_id(id_to_index,id)
        index = find(id_to_index==id);
        if size(index,2)==0
                index=-1
        endif;
endfunction;

# compute two lines intersection in 2d space given 4 points
# assume the two lines are never parallel
# inputs:
#   p_a,p_aa : start and end points on line a each [x,y]
#   p_b,p_bb : start and end points on line b each [x,y]
# output:
#   p : point of intersection [x,y]

function p = get2dLinesIntersection(p_a,p_aa,p_b,p_bb)
        m_a = (p_aa(2)-p_a(2))/(p_aa(1)-p_a(1));
        m_b = (p_bb(2)-p_b(2))/(p_bb(1)-p_b(1));
        b_a = p_a(2) - m_a * p_a(1);
        b_b = p_b(2) - m_b * p_b(1);
        x = (b_b - b_a)/(m_a - m_b);
        y = m_a * x + b_a;
        p=[x y]';
endfunction;

function v_idx=poseMatrixIndex(pose_index, num_poses, num_landmarks)
        global pose_dim;
        global landmark_dim;
        if (pose_index>num_poses)
                v_idx=-1;
                return;
        endif;
        v_idx=1+(pose_index-1)*pose_dim;
endfunction;

function v_idx=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks)
        global pose_dim;
        global landmark_dim;
        if (landmark_index>num_landmarks)
                v_idx=-1;
                return;
        endif;
        v_idx=1 + (num_poses)*pose_dim + (landmark_index-1) * landmark_dim;
endfunction;

function R=rotation2D(theta)
        s=sin(theta);
        c=cos(theta);
        R=[c -s;
        s  c];
endfunction;

function T=v2t(v)
        T=eye(3);
        T(1:2,1:2)=rotation2D(v(3));
        T(1:2,3)=v(1:2);
endfunction;

function v=t2v(A)
        v(1:2, 1)=A(1:2,3);
        v(3,1)=atan2(A(2,1),A(1,1));
endfunction;

function v=flattenIsometryByColumns(T)
v=zeros(6,1);
v(1:4)=reshape(T(1:2,1:2),4,1);
v(5:6)=T(1:2,3);
endfunction

function T=unflattenIsometryByColumns(v)
  T=eye(3);
  T(1:2,1:2)=reshape(v(1:4),2,2);
  T(1:2,3)=v(5:6);
endfunction

function J=J_atan2(p)
        J = 1./(p(1:2)'*p(1:2)) * [-p(2) p(1)];
endfunction;

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
