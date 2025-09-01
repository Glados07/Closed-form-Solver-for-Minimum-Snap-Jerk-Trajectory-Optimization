// - This is a distilled version of trajectory optimizer
// - Depending on <s>, the optimizer can solve up to Minimum Snap (s = 4) 
//   optimization problem by closed-form method.
// - The optimizer does not contain collision checking
// - The optimizer does not contain time allocation. However, 
//   a simple velocity & acceleration limit based heuristic method is provided.
/* -------------------------------How to use:-------------------------------*/
// 1. Initialize the optimizer with <s>, i.e. the order of optimization.
//    (s = 4: Minimum Snap, s = 3: Minimum Jerk, or lower but not less than 1)
// 2. Set one-dimensional waypoints (std::vector<double>) by calling setWaypoints() function.
// 3. Set piece time (std::vector<double>) by calling setPieceTime() function.
// 4. Call traj_optimize() function to get the trajectory 
//  - the function will return polynomial params (Eigen::VectorXd) of size M*(D+1), 
//    where M is the number of trajectory segments, 
//    D is the optimal polynomial degree (D = 2*s - 1)

# pragma once
#include <ros/ros.h>

# include <cmath>
# include <vector>
# include <eigen3/Eigen/Dense>

std::vector<double> simpleTimeAllocate(const std::vector<double>& waypoints, 
                                    double vel_limit, double acc_limit) {
    int M = waypoints.size() - 1;
    std::vector<double> piece_time(M, 0);
    std::vector<double> distance_ls(M, 0);
    for (int i = 0; i < M; i++){
        distance_ls[i] = std::abs(waypoints[i+1] - waypoints[i]);
    }
    double min_AccTime = vel_limit / acc_limit;
    double min_dist = 1.0/2 * acc_limit*min_AccTime*min_AccTime;
    for (int i = 0; i < M; i++){
        if (distance_ls[i] / 2.0 < min_dist){
            piece_time[i] = 2.0 * sqrt(distance_ls[i] / acc_limit);
        } else{
            piece_time[i] = 2*min_AccTime + (distance_ls[i] - 2*min_dist) / vel_limit;
        }
    }
    return piece_time;
}

template<int s>
class TrajOptimizer{
  private:
    std::vector<double> waypoints;
    int N; // number of waypoints
    int M; // number of segments
    int D; // polynomial degree
    bool waypoints_init_flag;

    std::vector<double> piece_time;
    bool piece_time_init_flag;

    Eigen::VectorXd res_params;
    std::vector<int> param_idx;   // use to locate the start index of each param in res_params / matrix from 0 -> M-1

    Eigen::MatrixXd QMat, AMat, JMat;

    Eigen::VectorXd vec_d, vec_d_F, vec_d_P;
    int dim_d_F, dim_d_P;

  public:
    TrajOptimizer() { 
        D = 2 * s - 1; 
        waypoints_init_flag = false;
        piece_time_init_flag = false;
    }

    bool setWaypoints(const std::vector<double>& wps) {
        waypoints.clear();
        waypoints = wps;

        N = waypoints.size();
        M = N - 1;

        piece_time.resize(M,0);
        waypoints_init_flag = true;

        return true;
    }

    bool setPieceTime(const std::vector<double>& pt) {
        if (pt.size() != (unsigned int)M) {
            ROS_ERROR("Optimizer: Error, piece_time size mismatch");
            return false;
        }
        if (!waypoints_init_flag) {
            ROS_ERROR("Optimizer: Error, waypoints not initialized");
            return false;
        }
        piece_time.clear();
        param_idx.clear();
        
        piece_time = pt;
        for (int i = 0; i < M; i++) {
            param_idx.push_back(i * (D + 1));
        }

        piece_time_init_flag = true;
        return true;
    }

    Eigen::VectorXd traj_optimize() {
        if (!waypoints_init_flag) {
            ROS_ERROR("Optimizer: Error, waypoints not initialized");
            return Eigen::VectorXd();
        } else if (!piece_time_init_flag) {
            ROS_ERROR("Optimizer: Error, piece_time not initialized");
            return Eigen::VectorXd();
        }
        QMat = Eigen::MatrixXd::Zero(M * (D + 1), M * (D + 1));
        AMat = Eigen::MatrixXd::Zero(M * (D + 1), M * (D + 1));
        vec_d = Eigen::VectorXd::Zero(M * (D + 1));

        dim_d_F = N + 2*(s - 1)+(M - 1) * (s + 1);
        dim_d_P = (M - 1) * (s - 2);
        vec_d_F = Eigen::VectorXd::Zero(dim_d_F);
        vec_d_P = Eigen::VectorXd::Zero(dim_d_P);

        res_params = Eigen::VectorXd::Zero(M * (D + 1));
        
        QMatGeneration();
        AMatGeneration();

        JMat = AMat.inverse().transpose() * QMat * AMat.inverse();

        // Extract the top-right submatrix R_FP from JMat
        Eigen::MatrixXd R_FP = JMat.block(0, M * (D + 1) - (M - 1) * (s - 2),       // start row & col
                                          M * (s + 2) + s - 2, (M - 1) * (s - 2));  // size
        // Extract the bottom-right submatrix R_PP from JMat
        Eigen::MatrixXd R_PP = JMat.block(M * (D + 1) - (M - 1) * (s - 2), M * (D + 1) - (M - 1) * (s - 2), // start point
                                         (M - 1) * (s - 2), (M - 1) * (s - 2));                             // size

        vec_d_F = vec_d.head(M * (s + 2) + s - 2);
        vec_d_P = -R_PP.inverse() * R_FP.transpose() * vec_d_F;

        vec_d.tail((M - 1) * (s - 2)) = vec_d_P;
        res_params = AMat.inverse() * vec_d;

        // Finished optimization, clear flags for next iteration
        waypoints_init_flag = false;
        piece_time_init_flag = false;

        return res_params;
    }

  private:
    void QMatGeneration(){
        Eigen::MatrixXd Q_elem = Eigen::MatrixXd::Zero(D+1, D+1);
        for (int i = 0; i < M; i++){
            Q_elem = Eigen::MatrixXd::Zero(D+1, D+1);
            for(int r = 0; r < D+1; r++){
                for(int c = 0; c < D+1; c++) {
                    if(D - r - s < 0 || D - c - s < 0) continue;
                    Q_elem(r, c) = factorial(D - r)/factorial(D - r - s)*
                                   factorial(D - c)/factorial(D - c - s)*
                                   1.0 / (2*D - 2*s - r - c + 1)*pow(piece_time[i], 2*D - 2*s - r - c);
                }
            }
            QMat.block(i*(D+1), i*(D+1), D+1, D+1) = Q_elem;
        }
    }

    void AMatGeneration(){
        int row_idx = 0;

        // 1 - waypoint constraints
        // from row 0 - row M
        for (int i = 0; row_idx < M; row_idx++, i++) {   // M waypoints
            AMat.block(row_idx, param_idx[i], 1, D+1) = posT_vec(0).transpose();
            vec_d[row_idx] = waypoints[row_idx];
        }
        // last waypoint
        AMat.block(row_idx, param_idx[M-1], 1, D+1) = posT_vec(piece_time[M-1]).transpose();
        vec_d[row_idx] = waypoints[N-1];
        row_idx++;

        // 2 - start & end differential constraints of s-1 order
        //   set all to 0 (rest status)
        for(int i = 1; i <= s-1; i++){
            switch (i) {
            case 1:
                AMat.block(row_idx, param_idx[0], 1, D+1) = velT_vec(0).transpose();
                vec_d[row_idx] = 0;
                break;
            case 2:
                AMat.block(row_idx, param_idx[0], 1, D+1) = accT_vec(0).transpose();
                vec_d[row_idx] = 0;
                break;
            case 3:
                AMat.block(row_idx, param_idx[0], 1, D+1) = jerkT_vec(0).transpose();
                vec_d[row_idx] = 0;
                break;
            case 4:
                AMat.block(row_idx, param_idx[0], 1, D+1) = snapT_vec(0).transpose();
                vec_d[row_idx] = 0;
                break;
            default:
                break;
            }
            row_idx++;
        }
        for(int i = 1; i <= s-1; i++){
            switch (i) {
            case 1:
                AMat.block(row_idx, param_idx[M-1], 1, D+1) = velT_vec(piece_time[M-1]).transpose();
                vec_d[row_idx] = 0;
                break;
            case 2:
                AMat.block(row_idx, param_idx[M-1], 1, D+1) = accT_vec(piece_time[M-1]).transpose();
                vec_d[row_idx] = 0;
                break;
            case 3:
                AMat.block(row_idx, param_idx[M-1], 1, D+1) = jerkT_vec(piece_time[M-1]).transpose();
                vec_d[row_idx] = 0;
                break;
            case 4:
                AMat.block(row_idx, param_idx[M-1], 1, D+1) = snapT_vec(piece_time[M-1]).transpose();
                vec_d[row_idx] = 0;
                break;
            default:
                break;
            }
            row_idx++;
        }

        // 3 - Continuity constraints of middle waypoints from order 0->s
        for(int pt_idx = 1; pt_idx <= M-1; pt_idx++){
            for(int i = 0; i <= s; i++){
                switch (i) {
                    case 0:
                        AMat.block(row_idx, param_idx[pt_idx-1], 1, D+1) = posT_vec(piece_time[pt_idx-1]).transpose();
                        AMat.block(row_idx, param_idx[pt_idx], 1, D+1) = -posT_vec(0).transpose();
                        vec_d[row_idx] = 0;
                        break;
                    case 1:
                        AMat.block(row_idx, param_idx[pt_idx-1], 1, D+1) = velT_vec(piece_time[pt_idx-1]).transpose();
                        AMat.block(row_idx, param_idx[pt_idx], 1, D+1) = -velT_vec(0).transpose();
                        vec_d[row_idx] = 0;
                        break;
                    case 2:
                        AMat.block(row_idx, param_idx[pt_idx-1], 1, D+1) = accT_vec(piece_time[pt_idx-1]).transpose();
                        AMat.block(row_idx, param_idx[pt_idx], 1, D+1) = -accT_vec(0).transpose();
                        vec_d[row_idx] = 0;
                        break;
                    case 3:
                        AMat.block(row_idx, param_idx[pt_idx-1], 1, D+1) = jerkT_vec(piece_time[pt_idx-1]).transpose();
                        AMat.block(row_idx, param_idx[pt_idx], 1, D+1) = -jerkT_vec(0).transpose();
                        vec_d[row_idx] = 0;
                        break;
                    case 4:
                        AMat.block(row_idx, param_idx[pt_idx-1], 1, D+1) = snapT_vec(piece_time[pt_idx-1]).transpose();
                        AMat.block(row_idx, param_idx[pt_idx], 1, D+1) = -snapT_vec(0).transpose();
                        vec_d[row_idx] = 0;
                        break;
                    default: break;
                }
                row_idx++;
            }
        }
        // 4 - free variables of middle points from order 1 --> s - 2
        for(int pt_idx = 1; pt_idx <= M-1; pt_idx++){
            for(int i = 1; i <= s - 2; i++){
                switch (i) {
                    case 1:
                        AMat.block(row_idx, param_idx[pt_idx], 1, D+1) = velT_vec(0).transpose();
                        // There is no need to set vec_d[row_idx] for free variables
                        break;
                    case 2:
                        AMat.block(row_idx, param_idx[pt_idx], 1, D+1) = accT_vec(0).transpose();
                        break;
                    case 3:
                        AMat.block(row_idx, param_idx[pt_idx], 1, D+1) = jerkT_vec(0).transpose();
                        break;
                    case 4:
                        AMat.block(row_idx, param_idx[pt_idx], 1, D+1) = snapT_vec(0).transpose();
                        break;
                    default: break;
                }
                row_idx++;
            }
        }

    }
  private:  //Misc Functions
    int factorial(int n){
        if(n <= 1) return 1;
        return n * factorial(n-1);
    }

    Eigen::VectorXd posT_vec(const double& t){
        Eigen::VectorXd res = Eigen::VectorXd::Zero(D+1);
        for(int i = 0; i <= D; i++){
            res(i) = pow(t, D - i);
        }
        return res;
    }

    Eigen::VectorXd velT_vec(const double& t){
        Eigen::VectorXd res = Eigen::VectorXd::Zero(D + 1);
        for(int i = 0; i <= D; i++){
            if(D-i-1 < 0) break;
            res(i) = (D - i) * pow(t, D - i - 1);
        }
        return res;
    }

    Eigen::VectorXd accT_vec(const double& t){
        Eigen::VectorXd res = Eigen::VectorXd::Zero(D + 1);
        for(int i = 0; i <= D; i++){
            if(D-i-2 < 0) break;
            res(i) = (D - i) * (D - i - 1) * pow(t, D - i - 2);
        }
        return res;
    }

    Eigen::VectorXd jerkT_vec(const double& t){
        Eigen::VectorXd res = Eigen::VectorXd::Zero(D + 1);
        for(int i = 0; i <= D; i++){
            if(D-i-3 < 0) break;
            res(i) = (D - i) * (D - i - 1) * (D - i - 2) * pow(t, D - i - 3);
        }
        return res;
    }

    Eigen::VectorXd snapT_vec(const double& t){
        Eigen::VectorXd res = Eigen::VectorXd::Zero(D + 1);
        for(int i = 0; i <= D; i++){
            if(D-i-4 < 0) break;
            res(i) = (D - i) * (D - i - 1) * (D - i - 2) * (D - i - 3) * pow(t, D - i - 4);
        }
        return res;
    }

};