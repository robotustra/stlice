My attempt to make fast console slicer for STL files.

15.07.2014
----------

The final goal.

	1) Stlice should read ascii STL files.
	2) Stlice is a console application.
	3) STLice generate contours and final path to generate g-code.


Input: 
	STL file
	Print parameters


STL file format:

solid name

facet normal ni nj nk
    outer loop
        vertex v1x v1y v1z
        vertex v2x v2y v2z
        vertex v3x v3y v3z
    endloop
endfacet

endsolid name

The main algorythm:

	1) Load points + projection of normales (could be up to 3 normales for each point)
	1.5) Load connection between points. (Joints.)
	2) Sort all vortex by Z.
	3) Find min max for x,y,z for all vortex.
	
	Loop to create contours:
	4) Start from minimal Z.
	5) Look for all joints to the upper dots.
	6) Find intersection of all joints with horizontal plane.
	7) Save layer, update connection between existing points and newly created.
	8) Proceed for all layers.

If Z spacing between points is too high, reduce the mesh along Z axis before start of procedure
