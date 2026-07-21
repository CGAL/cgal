
We provide various examples from the paper:
* [**Mesh projection and smoothing**](#mesh-smoothing-and-projection) with Geogram or Cinolib
* [**Usage of various queries**](#various-queries) with Geogram 
* [**Handle based deformation**](#handle-based-deformation) with LibIGL
* [**Offset computation with Alpha-wrapper**](#mesh-offsetting) with CGAL

We will slowly expand this list to match all figures in the paper.  <!-- for the [Graphics Replicability Stamp Initiative](https://www.replicabilitystamp.org/#). -->


# Mesh smoothing and projection
<img src="/assets/images/fitting.jpg" alt="Fitting image" width="600"/>

This executable is the most standard use of the code: improve the quality of mesh and fit to its geometric target. The code can be called as:
```Bash
./bin/geogram_smooth ../data/sphere.mesh ../data/max-planck.obj
```
*sphere.mesh* represent the volumetric tetrahedral mesh and *max-planck.obj* the geometric target represented as a surface triangle mesh. The code works with mesh containing Tetrahedra, Pyramids, Wedges and Hexahedra.

The same executable can also be called with a unique input parameter, as:
```Bash
./bin/geogram_smooth ../data/fandisk_kenshi_hexmesh.mesh
```
In this case it will operate as a smoothing (or untangling algorithm), improving the mesh quality while trying to preserve the outer boundary.

In both cases, the result will be saved as `output.mesh`. The executable `cinolib_smooth` operates in the exact same way and only differentiates from `geogram_smooth` because it uses a different frontend for mesh handling.

# Various queries
<img src="/assets/images/constraints.jpg" alt="Constraints image" width="900"/>

This executable highlight the alignment to all type of queries at the same time: Surface, Curves and a handle. The code can be called as:
```Bash
./bin/geogram_smooth_custom ../data/cube.mesh ../data/cube.mesh ../data/cube_target_cuves.mesh ../data/cube_target_points.txt
```
*cube.mesh* represent the input mesh, from which we read the cells (volume), the facets (surface) and the edges (curves). The second mesh is the surface used as projection, for which we simply use the same model. The third input "cube_target_cuves.mesh" contains the edge mesh forming the triangle that is used for projection. And finally, "cube_target_points.txt" contains vertex handle as: `id x y z`.

The final result will be saved as `output.mesh`. As mentioned in the paper, for this specific case, we use $\lambda_\text{strong} = 10$, see

```C++
    smoother.set_boundary_weight(Mesh_smoothing_3::Parameters::STRONG);
```
Your application may need a different parameter. 

# Handle based deformation
<img src="/assets/images/bunny.jpg" alt="Bunny handles" width="1200"/>


This executable will perform a deformation of the input object to match a set of input constraints, as shown in Fig. 12.
```Bash
./bin/libigl_handle_deformation ../data/bunny.mesh ../data/bunny_handles.txt
```
*bunny.mesh* is a standard tetrahedral mesh (as supported by LibIGL) and *bunny_handles.txt* contains vertex handle constraints as: `id x y z`.

Result will be saved as `handle_output.mesh`.

# Mesh offsetting

This executable can be used to reproduce Feature preserving offsets presented in the paper (Fig 14. and Fig 15.) We rely on the [Alpha-wrapper](https://doc.cgal.org/latest/Alpha_wrap_3/index.html) proposed by the CGAL library and use their parameters as input:
```Bash
./bin/cgal_offset_computation ../data/stairs.off 20 600
```
The input mesh must in the [OFF format](https://en.wikipedia.org/wiki/OFF_(file_format)) due to CGAL requirements. X=20 is used to compute the gate size in Alpha-wrap as alpha=D/X with D the diagonal of the bbox of the object. Y=600 is used to compute the offset of the mesh with delta=D/Y.

The executable will output two meshes: `offset_wrapper.mesh` and `offset_ours.mesh` to compare the before/after.


The parameters to run the different examples are as follows:
* Fig 14:
```Bash
./bin/cgal_offset_computation ../data/stairs.off 20 600
./bin/cgal_offset_computation ../data/stairs.off 40 600
./bin/cgal_offset_computation ../data/stairs.off 60 600
./bin/cgal_offset_computation ../data/stairs.off 80 600
```
<img src="/assets/images/stairs.jpg" alt="Stairs" width="1200"/>

* Fig 15:
```Bash
./bin/cgal_offset_computation ../data/43149.off 150 100 #  1% offset - slow remeshing (larger offset means slower)
./bin/cgal_offset_computation ../data/43149.off 150 50  #  2% offset
./bin/cgal_offset_computation ../data/43149.off 150 20  #  5% offset
./bin/cgal_offset_computation ../data/43149.off 150 10  # 10% offset
```
<img src="/assets/images/various_offsets.jpg" alt="Offset" width="1200"/>

The main bottleneck of the code is the inside remeshing that we perform with CGAL after the Alpha-wrapper.

