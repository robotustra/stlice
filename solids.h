/*
* This file contains loaded models data
*/


double * s_point[3]; 	// array of surface coordinates of model, X,Y,Z.
long   * pt_bond[2]; 	// one point can belong to many triangles. but we don't care about triangles.	
					 	//We care only about 
					 	// bonds. The number of bonds is 2 times less than the number of points.