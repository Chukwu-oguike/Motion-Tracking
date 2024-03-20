kalmanFilter: this function is an implementation of a kalman filter used to fuse sensor data from an accelerometer and 
a magnetometer with similar data dynamically calculated from a pivot frame of reference using a rate gyro. 

QuestAlgorithm: this is a function that uses a sequence of vectors from two separate reference frames
to determine the rotation matrix between both reference frames. The rotation matrix is the solution of a convex 
optimization problem.

FirstModel: this function applies the kalman filter to sensor data and utilizes the QUEST 
algorithm to determine the rotation matrix and quaternions representing
the relative rotations between two reference frames. It also plots the 
trajectory of the body undergoing the rotation.