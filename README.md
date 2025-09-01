# Closed form Solver for Minimum-Snap/Jerk Trajectory Optimization
This is a distilled version of the original closed-form method introduced by Charles Richter, Adam Bry, and Nicholas Roy.
The original paper: [Polynomial Trajectory Planning for Aggressive Quadrotor Flight in Dense Indoor Environments](https://link.springer.com/chapter/10.1007/978-3-319-28872-7_37)

## How to use
0. Include the TrajOptimizer.hpp file.
1. Instantiate an optimizer by the order of optimization (up to s = 4 for minimum snap). For example:
```
TrajOptimizer<4> agent;
```
2. Set **one-dimensional** waypoints (std::vector\<double\>) and segments time(std::vector\<double\>). Waypoints count N and segments time count M should satisfy M = N - 1. For example:
```
\\ wps has N elements
\\ piece_time has M = N - 1 elements
agent.setWaypoints(wps);
agent.setPieceTime(piece_time);
```
3. Call traj_optimize(). The agent will return the optimal polynomial coefficients (Eigen::VectorXd) with M*(D+1) elements (D is the optimal order of polynomial (D = 2*s-1)). For example:
```
Eigen::VectorXd poly_params;
poly_params = agent.traj_optimize();
```
![display](https://github.com/Glados07/Closed-form-Solver-for-Minimum-Snap-Jerk-Trajectory-Optimization/blob/main/display.gif)
## Acknowledgment
Sincere gratitude to repo [GCOPTER](https://github.com/ZJU-FAST-Lab/GCOPTER) which facilitated the verification and visualization of this code.
