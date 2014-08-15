/*
* This file contains loaded models data
*/


double * px; 	// array of surface coordinates of model, X,Y,Z.
double * py;
double * pz;
long   * pbond1; 	// one point can belong to many triangles. but we don't care about triangles.	
long   * pbond2; 	//We care only about 
				 	// bonds. The number of bonds is 2 times less than the number of points.
long * pt_face1;	// faces of points.
long * pt_face2;
long * pt_face3;

double * npx;	// normales.
double * npy;
double * npz;


long 	n_pt;	// the number of point in the model.
long 	n_uniq;
long 	n_pt1;	// the number of point to be added to the model after slicing objects.

// I think the same number of points should be generated for infill section on the base of normales information.
// It's kind of another nested solids with own boundaries.