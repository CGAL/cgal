# Migration from Surface_mesher to Mesh_3

Since `Surface_mesher` is deprecated, users should transition to `Mesh_3`. This package handles both surface and volume meshing.

Below is a comparison of the workflow for meshing an implicit function.

## Code Comparison

### 1. Headers and Types

**Old Workflow (`Surface_mesher`)**
```cpp
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
```

**New Workflow (`Mesh_3`)**
```cpp
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>

typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
```

### 2. Defining the Domain

**Old Workflow**
```cpp
// Defined by function pointer and bounding sphere
Surface_3 surface(sphere_function, Sphere_3(CGAL::ORIGIN, 2.)); 
```

**New Workflow**
```cpp
// Defined by a function wrapper (for labels) and bounding sphere
Mesh_domain domain(Function_wrapper(v), K::Sphere_3(CGAL::ORIGIN, 5.));
```

### 3. Meshing Criteria

**Old Workflow**
```cpp
CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30., 0.1, 0.1);
```

**New Workflow**
```cpp
// Criteria are split into Facet (surface) and Cell (volume)
// Args: angle, size, approximation
Facet_criteria facet_criteria(30, 0.2, 0.02); 
Mesh_criteria criteria(facet_criteria); 
```

### 4. Generation Function

**Old Workflow**
```cpp
CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
```

**New Workflow**
```cpp
C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, params::no_perturb());
```
## Common Pitfalls

* **Domain Wrapper:** Mesh_3 requires wrapping implicit functions (even single ones) to handle sub-domain labels correctly.

* **Squared Radius:** careful when porting bounding spheres, Mesh_3 examples often use K::Sphere_3 with a squared radius, whereas Surface_mesher examples might differ depending on the kernel.
