# eLife_2020_Analysis
Analysis code for calcium imaging data from Jain, Murphy-Baum, et al. (2020) eLife. 

This code contains several main functions, and supporting functions for use with analyzing calcium imaging data. These functions rely on a custom user interface, so modification will be necessary in order for the function to receive the right input parameters. However, the comments should make it clear enough what the intent of each code section is designed to do. 

Main functions: 

NT_GetROI() - extracts time-varying calcium imaging data from pre-defined ROIs.

AT_Correlate() - correlates one set of traces with another set, outputs the cross-correlograms.

VectorSum() - returns the vector sum angle or direction selective index from input tuning curves.

AT_AngularStats() - returns angular statistical information about the input traces. 
