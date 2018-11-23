# Clustering of non stationary time series with quadratic programming
## Abstract
The extraction of meaningful features from noised signals to characterize a system is one of the applications of time series analysis. Many techniques may be used for time series data mining. A possible approach to this problem is based on Clustering on Finite Element Method. From the point of view of practical implementation, the key ingredient of this method is the solution of a large-scale minimization problem.
This optimization problem is carried out iteratively by solving two consecutive convex subproblems in each step. The aim of this project is to investigate the performance of an Interior Point Method as solver for the computationally most expensive quadratic programming subproblem.
The algorithm performance is discussed by solving noised time series clustering problems.
## How to run
First, download the whole folder and open it in Matlab (Add to path -> selected folders and subfolders)

You can run some tests to see our implementation in action. Be aware that some tests might take a long time to finish.

The first test you can run is test_ipm.m, which is the core of the project. test_ipm.m tries to recover from a noised input signal the original data.

Other tests are test_time.m, which I used to test the performance, test_epsilon.m, which I used to observe the relation between the regularization parameter and the norm of the absolute error.
Also I included all the tests that I used for the plots in the report.
