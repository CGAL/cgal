#include <CGAL/Cartesian.h>
#include <CGAL/Conic_2.h>
#include <Moebius_diagram_2.h>
#include <Moebius_diagram_euclidean_traits_2.h>

typedef double NT;
typedef CGAL::Cartesian<NT> K;

struct Point : public K::Point_2; 

typedef double W;
typedef CGAL::Moebius_diagram_euclidean_traits_2<K,W> Traits;

struct MD : public CGAL::Moebius_diagram_2<Traits>;
struct WPoint : public MD::Point;
struct Vertex_iterator : public MD::Vertex_iterator;
struct Finite_edges_iterator : public MD::Finite_edges_iterator;
