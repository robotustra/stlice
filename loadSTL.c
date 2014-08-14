/*
 * 
 */

#include <string.h>
#include <stdio.h>

#include "loadSTL.h"


/* Custom fgets function to support Mac line-ends in Ascii STL files. OpenSCAD produces this when used on Mac */
void* fgets_(char* ptr, size_t len, FILE* f)
{
    while(len && fread(ptr, 1, 1, f) > 0)
    {
        if (*ptr == '\n' || *ptr == '\r')
        {
            *ptr = '\0';
            return ptr;
        }
        ptr++;
        len--;
    }
    return NULL;
}

/*
* 	1) First read just count the number of verteces in the file.
*	2) Allocate space for arrays.
*/

// This function returns the number of points in the model.
long countPointsSTLModel_ascii(const char* filename)
{
    FILE* f = fopen(filename, "rt");
    if (f == NULL) {
    	printf("No such file.\n");
    	return 0; // no such file
    }
    char buffer[1024];
    char pc[5] = "-/|\\";
    int n = 0;
    long nt = 0;
    double x,y,z;
    while(fgets_(buffer, sizeof(buffer), f))
    {
        if (sscanf(buffer, " vertex %lf %lf %lf", &x, &y, &z) == 3)
        {
        	printf("\rLoading ... %c", pc[nt%4]);
            n++;
            switch(n)
            {
            case 1:
                break;
            case 2:
                break;
            case 3:
                //vol->addFace(v0, v1, v2);
                n = 0;
                break;
            }
            nt++;
        }
    }
    printf("\nDone\n");
    fclose(f);
    return nt;
}


long loadModelSTL_ascii(const char* filename, long n_pt, double * px, double * py, double * pz,
	long * pt_face1, long * pt_face2, long * pt_face3, long * n_faces, long * n_faces_unique )
{
    FILE* f = fopen(filename, "rt");
    char buffer[1024];
    long ipt = 0, nt = 0;
    long i, j;
	double x,y,z;    
	int n = 0;
	long pt_index = 0;
	long ipf = 0; // count faces
	long ipf_uniq = 0; // unique faces
	long v1 = 0, v2 = 0, v3 = 0;

    while(fgets_(buffer, sizeof(buffer), f))
    {
        if (sscanf(buffer, " vertex %lf %lf %lf", &x, &y, &z) == 3)
        {
        	nt++;
        	// check unique points
        	j=0; // number of repetition
        	printf("\rAnalysing: %d%%", (int) ((nt*100)/n_pt) );
        	
            for (i = 0; i<ipt; i++)
            {
            	if ( (px[i] == x) && (py[i] == y) && (pz[i] == z ) ){
            		j++;
            		break; // i contains the point index.
            	}
            }
            if (j == 0) 
            { 	
            	px[ipt] = x;
	        	py[ipt] = y;
    	    	pz[ipt] = z;
    	    	ipt ++;
    	    	pt_index = ipt - 1;
    	    }
    	    else
    	    {
    	    	pt_index = i;
    	    }

           	// faces
        	n++;
            switch(n)
            {
            case 1:
            	v1 = pt_index;
                break;
            case 2:
            	v2 = pt_index;
                break;
            case 3:
            	v3 = pt_index;
            	ipf++;
            	//adding face if unique
                {
                	j=0;
                	// count unique faces 
		            for (i = 0; i<ipf_uniq; i++)
		            {
		            	if ( 	((pt_face1[i] == v1) && (pt_face2[i] == v2) && (pt_face3[i] == v3 )) ||
		            	 		((pt_face1[i] == v2) && (pt_face2[i] == v3) && (pt_face3[i] == v1 )) ||
		            	 		((pt_face1[i] == v3) && (pt_face2[i] == v1) && (pt_face3[i] == v2 )) ||
		            	 		((pt_face1[i] == v1) && (pt_face2[i] == v3) && (pt_face3[i] == v2 )) ||
		            	 		((pt_face1[i] == v3) && (pt_face2[i] == v2) && (pt_face3[i] == v1 )) ||
		            	 		((pt_face1[i] == v2) && (pt_face2[i] == v1) && (pt_face3[i] == v3 )) ){
		            		j++;
		            		break; // i contains the point index.
		            	}
		            }
		            if (j == 0) 
		            { 	
		            	pt_face1[ipf_uniq] = v1;
			        	pt_face2[ipf_uniq] = v2;
		    	    	pt_face3[ipf_uniq] = v3;
		    	    	ipf_uniq ++;
		    	    }
                }

                n = 0;
                break;
            }
        }
    }
    printf("\nDone\n");
    *n_faces = ipf;
    *n_faces_unique = ipf_uniq;
    fclose(f);
    return ipt; // unique number of points
}

/*
SimpleModel* loadModelSTL_binary(SimpleModel *m,const char* filename, FMatrix3x3& matrix)
{
    FILE* f = fopen(filename, "rb");
    char buffer[80];
    uint32_t faceCount;
    //Skip the header
    if (fread(buffer, 80, 1, f) != 1)
    {
        fclose(f);
        return nullptr;
    }
    //Read the face count
    if (fread(&faceCount, sizeof(uint32_t), 1, f) != 1)
    {
        fclose(f);
        return nullptr;
    }
    //For each face read:
    //float(x,y,z) = normal, float(X,Y,Z)*3 = vertexes, uint16_t = flags
    m->volumes.push_back(SimpleVolume());
    SimpleVolume* vol = &m->volumes[m->volumes.size()-1];
    if(vol == nullptr)
    {
        fclose(f);
        return nullptr;
    }

    for(unsigned int i=0;i<faceCount;i++)
    {
        if (fread(buffer, sizeof(float) * 3, 1, f) != 1)
        {
            fclose(f);
            return nullptr;
        }
        float v[9];
        if (fread(v, sizeof(float) * 9, 1, f) != 1)
        {
            fclose(f);
            return nullptr;
        }
        Point3 v0 = matrix.apply(FPoint3(v[0], v[1], v[2]));
        Point3 v1 = matrix.apply(FPoint3(v[3], v[4], v[5]));
        Point3 v2 = matrix.apply(FPoint3(v[6], v[7], v[8]));
        vol->addFace(v0, v1, v2);
        if (fread(buffer, sizeof(uint16_t), 1, f) != 1)
        {
            fclose(f);
            return nullptr;
        }
    }
    fclose(f);
    return m;
}

SimpleModel* loadModelSTL(SimpleModel *m,const char* filename, FMatrix3x3& matrix)
{
    FILE* f = fopen(filename, "r");
    char buffer[6];
    if (f == nullptr)
        return nullptr;

    if (fread(buffer, 5, 1, f) != 1)
    {
        fclose(f);
        return nullptr;
    }
    fclose(f);

    buffer[5] = '\0';
    if (stringcasecompare(buffer, "solid") == 0)
    {
        SimpleModel* asciiModel = loadModelSTL_ascii(m, filename, matrix);
        if (!asciiModel)
            return nullptr;

        // This logic is used to handle the case where the file starts with
        // "solid" but is a binary file.
        if (m->volumes[m->volumes.size()-1].faces.size() < 1)
        {
            m->volumes.erase(m->volumes.end() - 1);
            return loadModelSTL_binary(m, filename, matrix);
        }
        return asciiModel;
    }
    return loadModelSTL_binary(m, filename, matrix);
}

SimpleModel* loadModelFromFile(SimpleModel *m,const char* filename, FMatrix3x3& matrix)
{
    const char* ext = strrchr(filename, '.');
    if (ext && strcmp(ext, ".stl") == 0)
    {
        return loadModelSTL(m,filename, matrix);
    }
    if (filename[0] == '#' && binaryMeshBlob != nullptr)
    {
        while(*filename == '#')
        {
            filename++;

            m->volumes.push_back(SimpleVolume());
            SimpleVolume* vol = &m->volumes[m->volumes.size()-1];
            int32_t n, pNr = 0;
            if (fread(&n, 1, sizeof(int32_t), binaryMeshBlob) < 1)
                return nullptr;
            cura::log("Reading mesh from binary blob with %i vertexes\n", n);
            Point3 v[3];
            while(n)
            {
                float f[3];
                if (fread(f, 3, sizeof(float), binaryMeshBlob) < 1)
                    return nullptr;
                FPoint3 fp(f[0], f[1], f[2]);
                v[pNr++] = matrix.apply(fp);
                if (pNr == 3)
                {
                    vol->addFace(v[0], v[1], v[2]);
                    pNr = 0;
                }
                n--;
            }
        }
        return m;
    }
    return nullptr;
}
*/
