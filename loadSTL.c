/*
 * 
 */

#include <string.h>
#include <stdio.h>

#include "solids.h"
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
    return nullptr;
}

/*
* 	1) First read just count the number of verteces in the file.
*	2) Allocate space for arrays.
*/

SimpleModel* loadModelSTL_ascii(SimpleModel *m,const char* filename, FMatrix3x3& matrix)
{
    m->volumes.push_back(SimpleVolume());
    SimpleVolume* vol = &m->volumes[m->volumes.size()-1];
    FILE* f = fopen(filename, "rt");
    char buffer[1024];
    FPoint3 vertex;
    int n = 0;
    Point3 v0(0,0,0), v1(0,0,0), v2(0,0,0);
    while(fgets_(buffer, sizeof(buffer), f))
    {
        if (sscanf(buffer, " vertex %lf %lf %lf", &vertex.x, &vertex.y, &vertex.z) == 3)
        {
            n++;
            switch(n)
            {
            case 1:
                v0 = matrix.apply(vertex);
                break;
            case 2:
                v1 = matrix.apply(vertex);
                break;
            case 3:
                v2 = matrix.apply(vertex);
                vol->addFace(v0, v1, v2);
                n = 0;
                break;
            }
        }
    }
    fclose(f);
    return m;
}
