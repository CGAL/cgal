//CGAL headers
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h> 

typedef CGAL::Simple_cartesian<double> Rp;
typedef Rp::Point_2 Point;
typedef Rp::Circle_2 Circle;
typedef double W;
typedef CGAL::Regular_triangulation_euclidean_traits_2<Rp,W>  Gt;
// typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
// typedef CGAL::Regular_triangulation_face_base_2<> Fb;
// typedef CGAL::Triangulation_data_structure_2<Vb,Fb > Tds;
typedef CGAL::Regular_triangulation_2<Gt> Regular_triangulation; 
typedef Regular_triangulation::Vertex_iterator Vertex_iterator;
