For the effective utilization of the algorithm P-SmvMDE algorithm, all three scripts have to be located in the same Matlab folder.

These three scripts are:

(1) P-SmvMDE_NCDF.m
(2) P-SmvMDE_NCDF_ms.m
(3) P-SmvMDE.m

The (1) P_SmvMDE_NCDF script is the first script that defines the operation of P-SmvMDE and for single time scale operations could be used independently even though it is not recommended.

The (2) P_SmvMDE_NCDF_ms.m script is a follow up script that slightly modifies the operation described in (1) to allow for the implementation of additional time scale operations while maintaining the original NCDF mapping as intented by design.

The (3) P_SmvMDE.m script is the recommended script to be used for calling the function since it provides complete utility both for single scale and multiscale operation by calling script (1) for the first scale and script (2) for any follow up scales.
