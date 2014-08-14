#include <stdio.h>
#include <stdlib.h>
#include "loadSTL.h"
#include "solids.h"


int main (int argc, char * argv[])
{
	long alloc_total = 0;
	long n_faces, n_faces_uniq;
	double min_x, max_x;
	double min_y, max_y;
	double min_z, max_z;
	double h_layer = 0.2; //default layers thickness if not overwritten by input
	long i, j; // counter
	long * ptl_down; // the points coordinates in the bottom of slab.
	long * ptl_up; // the points to which bottom points are connected. (to find directing vector.)
	long npt_bot = 0; // the number of points on the bottom.
	long npt_up = 0; // the number of point which are higher.
	int n_layers = 0; // the number of layers.
	int current_layer = 0;
	double ** layer_x; // x coordinates of particular layer.
	double ** layer_y; // y coordinates of particular layer.
	double ** layer_z; // z coordinates of particular layer.
	long * layer_npt; // the number of points allocated in the layer.

	if (argc == 1) {
		printf("Usage: %s [OPTIONS] stl_file\n", argv[0]);
		return 0;
	}

	// total number of points described in the file.
	n_pt = countPointsSTLModel_ascii(argv[1], &min_x, &max_x, &min_y, &max_y, &min_z, &max_z); 
	printf("Points in the STL file: %ld\n", n_pt);
	printf("Z min value: %.3lf\n", min_z);
	printf("Object size: X: %.3lfmm, Y: %.3lfmm, Z: %.3lfmm\n", max_x - min_x, max_y - min_y, max_z - min_z);

	//Allocate all points of model.
	px = (double *) malloc(n_pt * sizeof(double)); 
	py = (double *) malloc(n_pt * sizeof(double));
	pz = (double *) malloc(n_pt * sizeof(double));

	alloc_total += 3 * n_pt * sizeof(double);
	
	pbond1 = (long *) malloc(n_pt * sizeof(long)); // the number of pairs, the number of connection is wrong.

	alloc_total += n_pt * sizeof(long);

	pt_face1 = (long *) malloc(n_pt/3 * sizeof(long));
	pt_face2 = (long *) malloc(n_pt/3 * sizeof(long));
	pt_face3 = (long *) malloc(n_pt/3 * sizeof(long));

	alloc_total += n_pt * sizeof(long);
	
	printf("Allocated: %ld bytes\n", alloc_total);

	// Calculate the number of layers
	n_layers = (max_z - min_z) / h_layer;
	n_layers++;
	printf("Number of layers: %d\n", n_layers );
	// Allocate pointers for countur points. They might have different coordinates than original one.
	layer_x = (double **) malloc( n_layers * sizeof(double*)) ;
	layer_y = (double **) malloc( n_layers * sizeof(double*)) ;
	layer_z = (double **) malloc( n_layers * sizeof(double*)) ;

	layer_npt = (long *) malloc( n_layers * sizeof(long)) ;
	
	alloc_total += 3 * n_layers * sizeof(double*) + n_layers * sizeof(long);
	
	printf("Allocated: %ld bytes\n", alloc_total);
	//After we will allocate all layer counturs and save points of to that arrays.


	//Load model into memory
	n_uniq = loadModelSTL_ascii(argv[1], n_pt, px, py, pz, pt_face1, pt_face2, pt_face3, &n_faces, &n_faces_uniq);
	printf("Unique points in the model: %ld\n", n_uniq);
	printf("Number of faces: %ld\n", n_faces);
	printf("Number of unique faces: %ld\n", n_faces_uniq);

	// Should already know min max values.
	// Shifting all points by Z_min to get the surface points.

	if (min_z != 0) {
		printf("Shifting model by Z min.\n");
	}

	
	//#2 The number of points we will use is n_uniq.
	// The number of faces which is used to look through faces ia n_faces_uniq.

	//Let's put slice thickness is 0.2 mm for tests, after it will be input from console.
	//Looking for all point in the slab (0 .. 0.2) 
	// actually we could qsort all points by Z to speed up by slicing, but this could be done after. 
	// now we will just look through all list of points.

	// First just estimate the number of point in the first layer. After allocate memory, to prevent waste of
	// memory.

	ptl_down = (long *) malloc( n_uniq * sizeof(long));
	alloc_total += n_uniq * sizeof(long);

	ptl_up = (long *) malloc( n_uniq * sizeof(long));
	alloc_total += n_uniq * sizeof(long);
	
	printf("Allocated: %ld bytes\n", alloc_total);

	npt_bot = 0;

	for (i=0; i<n_uniq; i++)
	{
		if(  (0 <= pz[i]) && (pz[i] <= h_layer) )
		{
			//our client
			ptl_down[npt_bot] = i; // the index of the point.
			npt_bot++;
		}
	}
	printf("The number of points in layer %ld: %ld\n", 1, npt_bot);

	// Looking for all neighbours which have Z more than current slab and which are connected with current bottom points.

	npt_up = 0;
	for (i=0; i<n_faces_uniq; i++)
	{
		for(j=0; j<npt_bot; j ++)
		{
			if ( pt_face1[i] == ptl_down[j]) {
				// checking for neighbours of this point, if they have bigger Z.
				if ( pz[pt_face2[i]] > h_layer ) npt_up ++; 
				if ( pz[pt_face3[i]] > h_layer ) npt_up ++;
			}
			if ( pt_face2[i] == ptl_down[j]) {
				// checking for neighbours of this point, if they have bigger Z.
				if ( pz[pt_face1[i]] > h_layer ) npt_up ++; 
				if ( pz[pt_face3[i]] > h_layer ) npt_up ++;
			}
			if ( pt_face3[i] == ptl_down[j]) {
				// checking for neighbours of this point, if they have bigger Z.
				if ( pz[pt_face1[i]] > h_layer ) npt_up ++; 
				if ( pz[pt_face2[i]] > h_layer ) npt_up ++;
			}
		}
	}
	printf("The number of points connected to layer %ld: %ld\n", 1, npt_up);

	//free all stuff

	free(pt_face1);
	free(pt_face2);
	free(pt_face3);
	free(pbond1);
	free(px);
	free(py);
	free(pz);

	return 0;
}