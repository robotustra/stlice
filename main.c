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
	long npt_bot = 0; // the number of points on the bottom.
	long npt_up = 0; // the number of point which are higher.
	int n_layers = 0; // the number of layers.
	long current_layer = 0;
	double current_z_bot = 0; // bottom and top Z of the layer.
	double current_z_top = 0;
	double ** layer_x; // x coordinates of particular layer.
	double ** layer_y; // y coordinates of particular layer.
	//double ** layer_z; // z coordinates of particular layer.
	double * cnt_x; // x coordinaes of points of one layer.
	double * cnt_y; // y coord.
	double * cnt_x_next; // x coordinaes of points of the next layer which appear because of intersection of segments with plain.
	double * cnt_y_next; // y coord.
	
	// We are not allocating Z because we calculate Z coordinates from the thickness of layer and the number of layer.
	long * layer_npt; // the number of points allocated in the layer.
	long * cnt_pt_index; // current layer index of parent points. to find faces in the next layer.

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

	// First just estimate the number of point in the first layer. After allocate memory, to prevent waste of
	// memory.
	

	// Loop over all layers.
	for (current_layer=0; current_layer < n_layers; current_layer++)
	{
		current_z_bot = current_layer * h_layer;
		current_z_top = current_z_bot + h_layer;

		npt_bot = 0;
		// The first time we are looking for points in bottom layer.
		if (current_layer == 0) 
		{
			for (i=0; i<n_uniq; i++)
			{
				if( (current_z_bot <= pz[i]) && (pz[i] < current_z_top) )
				{
					npt_bot++;
				}
			}
		} 
		else
		{
			// All other times we already have all points of the next layer in the
			// array cnt_pt_index[] (it's index of parent point) with coordinaxes cnt_x_next[]
			// and cnt_y_next[]. The number of points in the new layer is "npt_up".

			// We also have to replace cordinates of points if the point from array ptl_up[npt_up]
			// fits into the current layer.
			for (i=0; i<npt_up; i++)
			{
				if ( (current_z_bot <= pz[ptl_up[i]]) && (pz[ptl_up[i]] < current_z_top) )
				{
					//replace the index of point, because we are already reached the end of line segment
					// with the precision of layer thickness.
					cnt_pt_index[i] = ptl_up[i]; // index of upper point which fit in the layer.
					cnt_x_next[i] = px[ptl_up[i]];
					cnt_y_next[i] = py[ptl_up[i]];
				}
			}
			npt_bot = npt_up;
		}
		printf("The number of points in layer %ld: %ld\n", current_layer + 1, npt_bot);

		ptl_down = (long *) malloc( npt_bot * sizeof(long));

		alloc_total += npt_bot * sizeof(long);

		// Allocate X and Y coordinates of current countur for all points of the model.
		cnt_x = (double *) malloc (npt_bot * sizeof(double));
		cnt_y = (double *) malloc (npt_bot * sizeof(double));

		alloc_total += 2 * npt_bot * sizeof(double);

		printf("Allocated: %ld bytes\n", alloc_total);

		
		// make this arrays addressible
		layer_x[current_layer] = cnt_x;
		layer_y[current_layer] = cnt_y;
		layer_npt[current_layer] = npt_bot;
		// This information about the countur of the layer.
		// Now copy the information about all points.

		if (current_layer == 0)
		{
			npt_bot = 0; // reset counter
			for (i=0; i<n_uniq; i++)
			{
				if(  ( (current_z_bot - very_small_number ) <= pz[i]) && (pz[i] < current_z_top) )
				{
					//our client
					ptl_down[npt_bot] = i; // the index of the point.
					layer_x[current_layer][npt_bot] = px[i];
					layer_y[current_layer][npt_bot] = py[i];
					npt_bot++;
				}
			}
		}
		else
		{
			//the next time we are already know all points we just need to copy it from cnt_pt_index[] array.
			for (i=0; i<npt_bot; i++)
			{
				ptl_down[i] = cnt_pt_index[i]; // index of upper point which fit in the layer.
				layer_x[current_layer][i] = cnt_x_next[i];
				layer_y[current_layer][i] = cnt_y_next[i];
			}
			// now we have all point in the current layer and the number of point is npt_bot.
		}

		printf("The number of points in layer %ld: %ld\n", current_layer + 1, npt_bot);

		// Looking for all neighbours which have Z more than current slab and which are connected with current bottom points.
		// The number of connection with upper layer gives the points of the upper layer.
		npt_up = 0;
		for (i=0; i<n_faces_uniq; i++)
		{
			for(j=0; j<npt_bot; j ++)
			{
				if ( pt_face1[i] == ptl_down[j]) {
					// checking for neighbours of this point, if they have bigger Z.
					if ( pz[pt_face2[i]] > current_z_top ) npt_up++; 
					if ( pz[pt_face3[i]] > current_z_top ) npt_up++;
				}
				if ( pt_face2[i] == ptl_down[j]) {
					// checking for neighbours of this point, if they have bigger Z.
					if ( pz[pt_face1[i]] > current_z_top ) npt_up++; 
					if ( pz[pt_face3[i]] > current_z_top ) npt_up++;
				}
				if ( pt_face3[i] == ptl_down[j]) {
					// checking for neighbours of this point, if they have bigger Z.
					if ( pz[pt_face1[i]] > current_z_top ) npt_up++; 
					if ( pz[pt_face2[i]] > current_z_top ) npt_up++;
				}
			}
		}
		printf("The number of points connected to layer %ld: %ld\n", current_layer + 1, npt_up);


		ptl_up = (long *) malloc( npt_up * sizeof(long));
		alloc_total += npt_up * sizeof(long);

		// Allocate X and Y coordinates of the next layer which will be generated. We know
		// that the number of points in the next countour will be "npt_up".
		cnt_x_next = (double *) malloc (npt_up * sizeof(double));
		cnt_y_next = (double *) malloc (npt_up * sizeof(double));
		cnt_pt_index = (long *) malloc (npt_up * sizeof(long));

		alloc_total += 2 * npt_up * sizeof(double) + npt_up * sizeof(long);

		printf("Allocated: %ld bytes\n", alloc_total);

		// Loadin all indexii of  upper points
		npt_up = 0;
		for (i=0; i<n_faces_uniq; i++)
		{
			for(j=0; j<npt_bot; j ++)
			{
				if ( pt_face1[i] == ptl_down[j]) {
					// checking for neighbours of this point, if they have bigger Z.
					if ( pz[pt_face2[i]] > current_z_top ) 
					{ 	
						ptl_up[npt_up] = pt_face2[i]; 
						//calculate xy coordinates
						cnt_pt_index[npt_up] = ptl_down[j]; // the index of parent point.
						// coefficient of proportionality to calculate new points coordinates
						k = (current_z_top - pz[ptl_down[j]]) / (pz[ptl_up[npt_up]] - pz[ptl_down[j]]); // this difference should not be eqial 0, but after I'll add check just in case
						cnt_x_next[npt_up] = px[ptl_down[j]] + k * (px[ptl_up[npt_up]] - px[ptl_down[j]]);
						cnt_y_next[npt_up] = py[ptl_down[j]] + k * (py[ptl_up[npt_up]] - py[ptl_down[j]]);

						npt_up++; 
					} 
					if ( pz[pt_face3[i]] > current_z_top ) 
					{ 
						ptl_up[npt_up] = pt_face3[i]; 
						//calculate xy coordinates
						cnt_pt_index[npt_up] = ptl_down[j]; // the index of parent point.
						// coefficient of proportionality to calculate new points coordinates
						k = (current_z_top - pz[ptl_down[j]]) / (pz[ptl_up[npt_up]] - pz[ptl_down[j]]); // this difference should not be eqial 0, but after I'll add check just in case
						cnt_x_next[npt_up] = px[ptl_down[j]] + k * (px[ptl_up[npt_up]] - px[ptl_down[j]]);
						cnt_y_next[npt_up] = py[ptl_down[j]] + k * (py[ptl_up[npt_up]] - py[ptl_down[j]]);

						npt_up++; 
					}
				}
				if ( pt_face2[i] == ptl_down[j]) {
					// checking for neighbours of this point, if they have bigger Z.
					if ( pz[pt_face1[i]] > current_z_top ) 
					{ 
						ptl_up[npt_up] = pt_face1[i]; 
						//calculate xy coordinates
						cnt_pt_index[npt_up] = ptl_down[j]; // the index of parent point.
						// coefficient of proportionality to calculate new points coordinates
						k = (current_z_top - pz[ptl_down[j]]) / (pz[ptl_up[npt_up]] - pz[ptl_down[j]]); // this difference should not be eqial 0, but after I'll add check just in case
						cnt_x_next[npt_up] = px[ptl_down[j]] + k * (px[ptl_up[npt_up]] - px[ptl_down[j]]);
						cnt_y_next[npt_up] = py[ptl_down[j]] + k * (py[ptl_up[npt_up]] - py[ptl_down[j]]);

						npt_up++; 
					}
					if ( pz[pt_face3[i]] > current_z_top ) 
					{ 
						ptl_up[npt_up] = pt_face3[i]; 
						cnt_pt_index[npt_up] = ptl_down[j]; // the index of parent point.
						// coefficient of proportionality to calculate new points coordinates
						k = (current_z_top - pz[ptl_down[j]]) / (pz[ptl_up[npt_up]] - pz[ptl_down[j]]); // this difference should not be eqial 0, but after I'll add check just in case
						cnt_x_next[npt_up] = px[ptl_down[j]] + k * (px[ptl_up[npt_up]] - px[ptl_down[j]]);
						cnt_y_next[npt_up] = py[ptl_down[j]] + k * (py[ptl_up[npt_up]] - py[ptl_down[j]]);
					
						npt_up++; 
					}
				}
				if ( pt_face3[i] == ptl_down[j]) {
					// checking for neighbours of this point, if they have bigger Z.
					if ( pz[pt_face1[i]] > current_z_top ) 
					{ 
						ptl_up[npt_up] = pt_face1[i];
						cnt_pt_index[npt_up] = ptl_down[j]; // the index of parent point.
						// coefficient of proportionality to calculate new points coordinates
						k = (current_z_top - pz[ptl_down[j]]) / (pz[ptl_up[npt_up]] - pz[ptl_down[j]]); // this difference should not be eqial 0, but after I'll add check just in case
						cnt_x_next[npt_up] = px[ptl_down[j]] + k * (px[ptl_up[npt_up]] - px[ptl_down[j]]);
						cnt_y_next[npt_up] = py[ptl_down[j]] + k * (py[ptl_up[npt_up]] - py[ptl_down[j]]);

						npt_up++; 
					}
					if ( pz[pt_face2[i]] > current_z_top ) 
					{ 
						ptl_up[npt_up] = pt_face2[i]; 
						cnt_pt_index[npt_up] = ptl_down[j]; // the index of parent point.
						// coefficient of proportionality to calculate new points coordinates
						k = (current_z_top - pz[ptl_down[j]]) / (pz[ptl_up[npt_up]] - pz[ptl_down[j]]); // this difference should not be eqial 0, but after I'll add check just in case
						cnt_x_next[npt_up] = px[ptl_down[j]] + k * (px[ptl_up[npt_up]] - px[ptl_down[j]]);
						cnt_y_next[npt_up] = py[ptl_down[j]] + k * (py[ptl_up[npt_up]] - py[ptl_down[j]]);

						npt_up++; 
					}
				}
			}
		}

		// Now we know all points where bottom points go, and generated temporary points which intersect with
		// current_z_top plane.
		// We have save the point (x,y) coordinates and the "index" which comes from parent point.

	} //loop over all layers.

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