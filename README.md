# Presentation
This project aims to compare different methods in computing the trace of the inverse of a symmetric invertible matrix $A$ :<br />
<p align="center">
$$tr(A^{-1})$$
</p>
The main goal is to implement Algorithm 1 and 2 from [1], reproduce the results of [1], and compare the running times of Algorithm 1 and 2
with other basic algorithms to compute the trace of the inverse. This project follows the questions of Trace_estimator_questions.pdf.


# Explanation
The code of the project is organized as follows :

-```/src/adjust.m```, ```/src/Algorithm1.m``` and ```/src/Algorithm2.m``` are matlab function files used to implement Algorithm 1 and 2 from [1].<br />
-```/src/formVFH.m``` performs 1 step of the recurrence relation to construct the VFH matrix from [1].<br />
-```/src/Viscek_matrix.m```, ```/src/Heat_matrix.m``` and ```/src/Poisson_matrix.m``` produce the different tables from [1].<br />
-```/src/trace_invA.m``` and ```/src/time_measure.m``` produce Figure 1 and 2 respectively from [1]. Examples of the plots used to construct those figures are stored in ```/plots```.<br />

Note that the files used to produce the different results were executed on Matlab R2022b and do not need any special libraries.

# References
[1] Zhaojun Bai, Gark Fahey, and Gene Golub. Some large-scale matrix computation
problems. Journal of Computational and Applied Mathematics, 74(1):71–89, 1996.<br /><br />

[2] Gene H. Golub and John H. Welsch. Calculation of gauss quadrature rules. Mathematics
of Computation, 23(106):221–s10, 1969.
