
long countPointsSTLModel_ascii(const char* filename, 
	double * minx, double * maxx, double * miny, double * maxy, double * minz, double * maxz);
long loadModelSTL_ascii(const char* filename, long n_pt, double * px, double * py, double * pz,
	long * pt_face1, long * pt_face2, long * pt_face3, long * n_faces, long * n_faces_unique,
    long * pt_bond1, long * pt_bond2, long * n_bonds );