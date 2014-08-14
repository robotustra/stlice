#include <stdio.h>
#include <stdlib.h>
#include "loadSTL.h"
#include "solids.h"


int main (int argc, char * argv[])
{
	long alloc_total = 0;
	long n_faces, n_faces_uniq;

	if (argc == 1) {
		printf("Usage: %s [OPTIONS] stl_file\n", argv[0]);
		return 0;
	}

	// total number of points described in the file.
	n_pt = countPointsSTLModel_ascii(argv[1]); 
	printf("Points in the model: %ld\n", n_pt);

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

	//Load model into memory
	n_uniq = loadModelSTL_ascii(argv[1], n_pt, px, py, pz, pt_face1, pt_face2, pt_face3, &n_faces, &n_faces_uniq);
	printf("Unique points in the model: %ld\n", n_uniq);
	printf("Number of faces: %ld\n", n_faces);
	printf("Number of unique faces: %ld\n", n_faces_uniq);





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