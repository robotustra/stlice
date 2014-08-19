#include <stdio.h>
#include <stdlib.h>
#include "loadSTL.h"
#include "solids.h"
#include <math.h>

#define very_small_number	1.0e-10



int main (int argc, char * argv[])
{
	long alloc_total = 0;
	long n_faces, n_faces_uniq;
	double min_x, max_x;
	double min_y, max_y;
	double min_z, max_z;
	double h_layer = 0.2; //default layers thickness if not overwritten by input
	long i, j, l; // counter
	double k; // coefficient
	long * ptl_down; // the points coordinates in the bottom of slab.
	long * ptl_up; // the points to which bottom points are connected. (to find directing vector.)
	long npt_gen = 0; // the number of points being generated for layers.
	long npt_layer = 0; 
	long npt_total = 0;
	int n_layers = 0; // the number of layers.
	long current_layer = 0;
	long * current_layer_index;
	double current_z_bot = 0; // bottom and top Z of the layer.
	double current_z_top = 0;
	//double ** layer_x; // x coordinates of particular layer.
	//double ** layer_y; // y coordinates of particular layer.
	long ** layer_pt_index; // this is a 
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
	long * pt_bond1;
	long * pt_bond2;
	long n_bonds = 0;
	long n_start, n_finish;
	double x_start, x_finish;
	double y_start, y_finish;
	double z_start, z_finish;
	double ax,ay,az,bx,by,bz; // temporary vctors.	
	long norm_cnt = 0; // normales counter.
	long l_index = 0;
	double r_length = 0;
	long l_npt = 0; // the number of points in the layer.
	long p1_idx, p2_idx;



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
	
	//bonds
	pt_bond1 = (long *) malloc (n_pt / 2 * sizeof(long));
	pt_bond2 = (long *) malloc (n_pt / 2 * sizeof(long));

	alloc_total +=  n_pt * 0.5 * sizeof(long);	
	

	pt_face1 = (long *) malloc(n_pt/3 * sizeof(long));
	pt_face2 = (long *) malloc(n_pt/3 * sizeof(long));
	pt_face3 = (long *) malloc(n_pt/3 * sizeof(long));

	alloc_total += n_pt * sizeof(long);

	
	printf("Allocated: %ld bytes\n", alloc_total);

	// Calculate the number of layers
	n_layers = (max_z - min_z + very_small_number) / h_layer;
	n_layers++;
	printf("Number of layers: %d\n", n_layers );
	// Allocate pointers for countur points. They might have different coordinates than original one.
	//layer_x = (double **) malloc( n_layers * sizeof(double*)) ;
	//layer_y = (double **) malloc( n_layers * sizeof(double*)) ;
	layer_pt_index = (long ** ) malloc( n_layers * sizeof(long*)) ;

	//layer_z = (double **) malloc( n_layers * sizeof(double*)) ; // don't need Z coordinate because we evaluate Z from layer number.

	layer_npt = (long *) malloc( n_layers * sizeof(long)) ;
	
	alloc_total += n_layers * sizeof(long*) + n_layers * sizeof(long);
	
	printf("Allocated: %ld bytes\n", alloc_total);
	//After we will allocate all layer counturs and save points of to that arrays.

	//Load model into memory
	n_uniq = loadModelSTL_ascii(argv[1], n_pt, px, py, pz, pt_face1, pt_face2, pt_face3, &n_faces, &n_faces_uniq,
		pt_bond1, pt_bond2, &n_bonds);

	printf("Unique points in the model: %ld\n", n_uniq);
	printf("Number of faces: %ld\n", n_faces);
	printf("Number of unique faces: %ld\n", n_faces_uniq);
	printf("Number of unique pairs: %ld\n", n_bonds);

	// calculate normales to faces
	npx = (double *) malloc ( n_faces_uniq * sizeof(double));
	npy = (double *) malloc ( n_faces_uniq * sizeof(double));
	npz = (double *) malloc ( n_faces_uniq * sizeof(double));
	
	alloc_total += n_faces_uniq * 3 * sizeof(double);
	
	printf("Allocated: %ld bytes\n", alloc_total);

	// the number of normales is exactly the same as the number of faces.
	// And the index of normale coniside with the index of face.

	for (i=0; i<n_faces_uniq; i++)
	{
		ax = px[pt_face2[i]] - px[pt_face1[i]];
		ay = py[pt_face2[i]] - py[pt_face1[i]];
		az = pz[pt_face2[i]] - pz[pt_face1[i]];

		bx = px[pt_face3[i]] - px[pt_face2[i]];
		by = py[pt_face3[i]] - py[pt_face2[i]];
		bz = pz[pt_face3[i]] - pz[pt_face2[i]];

		npx[i] = ay*bz - az*by;
		npy[i] = az*bx - ax*bz;
		npz[i] = ax*by - ay*bx;
		r_length = sqrt( npx[i]*npx[i] + npy[i]*npy[i] + npz[i]*npz[i] );

		npx[i] /= r_length;
		npy[i] /= r_length;
		npz[i] /= r_length;

	}
	printf("Normales done\n");

	// Should already know min max values.
	// Shifting all points by Z_min to get the surface points.

#if 0
	//count how many normales have every point,
	for (i=0; i<n_uniq; i++)
	{
		norm_cnt = 0;
		for(j=0; j<n_faces_uniq; j++){
			if( (pt_face1[j] == i) || (pt_face2[j] == i) || (pt_face1[j] == i)) {
				norm_cnt ++;
			}
		}
		printf("Point %ld, normales: %ld\n", i, norm_cnt);
	}

	if (min_z != 0) {
		printf("Shifting model by Z min.\n");
	}
#endif
	
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
	/*
	for (i=0; i< n_bonds; i++)
	{
		//consider all points which has different Z, crossing the horizontal plane.
		v1 = pt_bond1[i];
		v2 = pt_bond2[i];
		npt_gen += (long)  ( (abs( pz[v1] - pz[v2] ) + very_small_number) / h_layer + 2 ); // two extra point to include ends,
	}*/
	npt_gen = 0;
	for (i=0; i< n_bonds; i++)
	{
		//consider all points which has different Z, crossing the horizontal plane.
		v1 = pt_bond1[i];
		v2 = pt_bond2[i];
		
		nl1 = pz[v1] / h_layer; // the number of layer where the point resides
		nl2 = pz[v2] / h_layer;

		if ( (nl1 - nl2) != 0 )
		{
			//have some points to add.
			if ( nl1 > nl2 ) 
			{ 
				n_start = nl2; 
				n_finish = nl1;
			}
			else
			{
				n_start = nl1;
				n_finish = nl2;
			}
			for (j = n_start+1; j< n_finish ; j++ ) npt_gen++;	
		}
	}
	printf("Number of points to generate: %ld\n", npt_gen );
	// adding points from initial array.

	// Now I have to generate all point and after distribute it over layers
	// The total amount of point

	npt_total = npt_gen + n_uniq; // it will be all points in the countur.

	cnt_x = (double *) malloc (npt_total * sizeof(double));
	cnt_y = (double *) malloc (npt_total * sizeof(double));
	cnt_z = (double *) malloc (npt_total * sizeof(double));

	alloc_total += 3 * npt_total * sizeof(double);
	
	printf("Allocated: %ld bytes\n", alloc_total);

	
	// Copy unique point to the common array for simplicity
	npt_gen = 0;
	for (i=0; i<n_uniq; i++){
		cnt_z[npt_gen] = pz[i];
		cnt_x[npt_gen] = px[i];
		cnt_y[npt_gen] = py[i];
		npt_gen++;
	}

#if 0	
	// and Generate all coordinates of our new points
	for (i=0; i< n_bonds; i++)
	{
		//consider all points which has different Z, crossing the horizontal plane.
		v1 = pt_bond1[i];
		v2 = pt_bond2[i];
		
		nl1 = pz[v1] / h_layer; // the number of layer where the point resides
		nl2 = pz[v2] / h_layer;

		if ( (nl1 - nl2) != 0 )
		{
			//have some points to add.
			// All Z coordinates of newly generated points will be maped on Z = layer_numb * h_layer; layer_numb: (0 .. n_layers)
			// Just calculate points coordinates form linear interpolation.
			if ( nl1 > nl2 ) 
			{ 
				n_start = nl2; 
				n_finish = nl1;
				x_start = px[v2];
				x_finish = px[v1];
				y_start = py[v2];
				y_finish = py[v1];
				z_start = pz[v2];
				z_finish = pz[v1];
			}
			else
			{
				n_start = nl1;
				n_finish = nl2;
				x_start = px[v1];
				x_finish = px[v2];
				y_start = py[v1];
				y_finish = py[v2];
				z_start = pz[v1];
				z_finish = pz[v2];
			}
			
			for (j = n_start+1; j< n_finish ; j++ )
			{
				cnt_z[npt_gen] = j * h_layer;
				k = (j * h_layer - z_start) / (z_finish - z_start); 
				cnt_x[npt_gen] = x_start + k * (x_finish - x_start);
				cnt_y[npt_gen] = y_start + k * (y_finish - y_start);

				npt_gen++;	
			}
		}
	}
#endif
	// Now we have all point in one array but the indexing from faces to points is saved, because

	// this is a real number of allocated poins
	printf("All points: %ld\n", npt_gen );
	
	// now we have all points, can distribute them in layers.
	npt_total = 0;
	for (j=0; j<n_layers; j++)
	{
		current_z_bot = j * h_layer;
		current_z_top = current_z_bot + h_layer;

		npt_layer = 0;
		// generated points
		for (i=0; i<npt_gen; i++)
		{
			if ( ((current_z_bot - very_small_number) <= cnt_z[i]) && 
				 ((current_z_top - very_small_number) > cnt_z[i])  )
			{
				npt_layer ++;
			}
		}
		// allocate counturs.

		current_layer_index = (long *) malloc (npt_layer * sizeof(long));

		npt_layer = 0;
		for (i=0; i<npt_gen; i++)
		{
			if ( ((current_z_bot - very_small_number) <= cnt_z[i]) && 
				 ((current_z_top - very_small_number) > cnt_z[i])  )
			{
				current_layer_index[npt_layer] = i;
				npt_layer ++;
			}
		}

		npt_total += npt_layer;
		layer_npt[j] = npt_layer;
		layer_pt_index[j] = current_layer_index;

#if 0
		printf("Layer: %ld, Points: %ld\n", j, npt_layer );
#endif

	}
	printf("Total points after layering: %ld\n", npt_total );

	//#3 Ready to process counturs.

#if 0
	//Test print of countur point.
	//output all points of the first layer.
	l_index = 0;
	for (i=0; i< layer_npt[l_index]; i++)
	{
		printf(" %ld %.3lf %.3lf %.3lf\n", i, cnt_x[layer_pt_index[l_index][i]],
			cnt_y[layer_pt_index[l_index][i]], cnt_z[layer_pt_index[l_index][i]] );
	}
#endif

	//ordering points accordint to normals information.

	l_index = 0;
	l_npt = layer_npt[l_index]; // the number of points in the layer.
	long ii = 0;
	for (i=0; i< l_npt; i++)
	{
		p1_idx = layer_pt_index[l_index][i];

		for (l=i+1; l<l_npt; l++){

			p2_idx = layer_pt_index[l_index][l];

			//x = cnt_x[p_idx];
			//y = cnt_y[p_idx];
			//z = cnt_z[p_idx];

			for (j=0; j<n_faces_uniq; j++)
			{
				if ( (pt_face1[j] == p1_idx && pt_face2[j] == p2_idx && pz[pt_face3[j]] > pz[p1_idx]) ||
					 (pt_face2[j] == p1_idx && pt_face3[j] == p2_idx && pz[pt_face1[j]] > pz[p1_idx]) ||
					 (pt_face3[j] == p1_idx && pt_face1[j] == p2_idx && pz[pt_face2[j]] > pz[p1_idx]) ||
					 (pt_face1[j] == p1_idx && pt_face3[j] == p2_idx && pz[pt_face2[j]] > pz[p1_idx]) ||
					 (pt_face3[j] == p1_idx && pt_face2[j] == p2_idx && pz[pt_face1[j]] > pz[p1_idx]) ||
					 (pt_face2[j] == p1_idx && pt_face1[j] == p2_idx && pz[pt_face3[j]] > pz[p1_idx]) 
					 ) 
				{
					// the bond in the plane is found.
					// it might happen that one vertex has more than two bonds in the 
					// plane Z=0.

					printf ("i: %ld, j: %ld\n", p1_idx, p2_idx);
					ii++;
				}
			}
		}
		//looking for this point in the list of faces

	}
	printf("Pints: %ld, Bonds: %ld\n", l_npt, ii);


	// Now I have to invent how to create a contour on points of layer



	//free all stuff

	free(pt_face1);
	free(pt_face2);
	free(pt_face3);
	free(pt_bond1);
	free(pt_bond2);
	free(px);
	free(py);
	free(pz);

	return 0;
}