#include <CGAL/AABB_intersections.h>
#include "Polyhedron_type.h"
#include "Image_type.h"
#include <CGAL/AABB_tree.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Labeled_image_mesh_domain_3.h>

typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, Kernel> Mesh_domain;
typedef CGAL::Labeled_image_mesh_domain_3<Image,Kernel> Image_mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;

// 3D complex
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
