#include <stdio.h>
#include <stdlib.h>
#include "loadSTL.h"
#include "solids.h"

#define very_small_number	1.0e-10



int main (int argc, char * argv[])
{
	long alloc_total = 0;
	long n_faces, n_faces_uniq;
	double min_x, max_x;
	double min_y, max_y;
	double min_z, max_z;
	double h_layer = 0.2; //default layers thickness if not overwritten by input
	long i, j; // counter
	double k; // coefficient
	long * ptl_down; // the points coordinates in the bottom of slab.
	long * ptl_up; // the points to which bottom points are connected. (to find directing vector.)
	long npt_gen = 0; // the number of points being generated for layers.
	long npt_total = 0;
	int n_layers = 0; // the number of layers.
	long current_layer = 0;
	double current_z_bot = 0; // bottom and top Z of the layer.
	double current_z_top = 0;
	double ** layer_x; // x coordinates of particular layer.
	double ** layer_y; // y coordinates of particular layer.
	//double ** layer_z; // z coordinates of particular layer.
	double * cnt_x; // x coordinaes of points of one layer.
	double * cnt_y; // y coord.
	double * cnt_z; // z coord.

	double * cnt_x_next; // x coordinaes of points of the next layer which appear because of intersection of segments with plain.
	double * cnt_y_next; // y coord.
	
	// We are not allocating Z because we calculate Z coordinates from the thickness of layer and the number of layer.
	long * layer_npt; // the number of points allocated in the layer.
	long * cnt_pt_index; // current layer index of parent points. to find faces in the next layer.
	long v1, v2, v3;
	long nl1, nl2, nl3;

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
	//layer_z = (double **) malloc( n_layers * sizeof(double*)) ; // don't need Z coordinate because we evaluate Z from layer number.

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

	// First just estimate the number of point which we create in all layers. After allocate memory, to prevent waste of
	// memory.

	// Do loop over all faces, and consider all points crossing section planes.
	npt_gen = 0;
	for (i=0; i< n_faces_uniq; i++)
	{
		//consider all points which has different Z, crossing the horizontal plane.
		v1 = pt_face1[i];
		v2 = pt_face2[i];
		v3 = pt_face3[i];
		npt_gen += (long) abs( pz[v1] - pz[v2] ) / h_layer;
		npt_gen += (long) abs( pz[v2] - pz[v3] ) / h_layer;
		npt_gen += (long) abs( pz[v1] - pz[v3] ) / h_layer;
	}

	npt_gen >>= 1; // divide by 2 because the pair of verteces are taken into account 2 times.
	
	printf("Number of points to generate: %ld\n", npt_gen );

	// Now I have to generate all point and after distribute it over layers
	// The total amount of point

	npt_total = n_uniq + npt_gen; // it will be all points in the countur.

	cnt_x = (double *) malloc (npt_total * sizeof(double));
	cnt_y = (double *) malloc (npt_total * sizeof(double));
	cnt_z = (double *) malloc (npt_total * sizeof(double));

	alloc_total += 3 * npt_total * sizeof(double);
	
	printf("Allocated: %ld bytes\n", alloc_total);

	//Generate all coordinates of our new points
	npt_gen = 0;
	for (i=0; i< n_faces_uniq; i++)
	{
		//consider all points which has different Z, crossing the horizontal plane.
		v1 = pt_face1[i];
		v2 = pt_face2[i];
		v3 = pt_face3[i];
		
		nl1 = pz[v1] / h_layer; // the number of layer where the point resides
		nl2 = pz[v2] / h_layer;
		nl3 = pz[v3] / h_layer;

		if ( (nl1 - nl2) != 0 )
		{
			//have some points to add.
			// All Z coordinates of newly generated points will be maped on Z = layer_numb * h_layer; layer_numb: (0 .. n_layers)
			// Just calculate points coordinates form linear interpolation.
			n_start = nl1;
			n_finish = nl2;
			if ( nl1 > nl2 ) { n_start = nl2; n_finish = nl1 };
			for (j = 0; j< abs(nl1 - nl2); j++ )
			{
				cnt_z[npt_gen] = (n_start + j) * h_layer;
				cnt_x[npt_gen] = ....
				cnt_y[npt_gen] = ....

				npt_gen++;	
			}
			  
		}
		
	}

	





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