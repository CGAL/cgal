Benchmarks to compare Linear_cell_complex with other libraries:

In 2D:

* Polyhedron_3 (from CGAL)
* OpenMesh
* SurfaceMesh
* CGoGN

In 3D:

* OpenVolumeMesh
* CGoGN

INSTALLATION:
==============

1) Install all the following libraries:

CGAL: https://www.cgal.org/
CGoGN: https://cgogn.github.io/
OpenMesh: https://www.openmesh.org/
OpenVolumeMesh: https://www.openvolumemesh.org/

2) create links (or copy directory):
* in the 2D directory:
    one link called cgogn through the install directory of CGoGN;
    one link called openmesh through the install directory of OpenMesh.

* in the 3D directory:
    one link called cgogn through the install directory of CGoGN;
    one link called openvolumemesh through the openvolumemesh directory itself.

3) Compile the benchmarks

* create a build directory; cd build; cmake -DCGAL_DIR=CGAL_BUILD_DIR ..

CGAL_BUILD_DIR being the build directory of the CGAL library.

4) run the benchmarks.
* In 2D, the programs take off files as input.
* In 3D, lcc and cgogn take tetmesh and OpenVolumeMesh takes ovm.

You can create a tetmesh file using tetgen with an OFF file as input (https://www.berlios.de/software/tetgen/) with option -g to generate XXX.mesh file. Rename this file into XXX.tetmesh. Modify the file to keep only the two following sections:

**********************
Vertices
26002
...
Tetrahedra
174445
...
**********************

And for each tetrahedra line, remove the last column.

You obtain the tetmesh file. Then you can use the FileConverter/file_converter program in the OpenVolumeMesh library to convert it into the corresponding ovm file.
