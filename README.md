My attempt to make fast console slicer for STL files.

15.07.2014
----------

The final goal.

	1) STLice should read ascii STL files.
	2) STLice is a console application.
	3) STLice generate contours and final path to generate g-code.
	4) STLice should scale STL model with independent scale factors along x,y,z axii.
	5) STLice should be able to slice multiple models to be printed at the same time.
		Each model could have own scale factor(s).
	6) Optional: Generate code for multiple extuders on the base of color information of model. 


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

The main algorithm:

	1) Load points + projection of normals (could be up to 3 normals for each point)
	1.5) Load connection between points. (Joints.)
	1.6) Scale models and place them accordingly if there are many models.
	2) Sort all vortex by Z.
	3) Find min max for x,y,z for all vortex.
	
	Loop to create contours:
	4) Start from minimal Z.
	5) Look for all joints to the upper dots.
	6) Find intersection of all joints with horizontal plane.
	7) Save layer, update connection between existing points and newly created.
	8) Proceed for all layers.

If Z spacing between points is too high, reduce the mesh along Z axis before start of procedure
