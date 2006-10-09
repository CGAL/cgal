#include <cstdio> 
#include <fstream>
#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

using namespace boost ;
using namespace std ;
using namespace CGAL ;

typedef Simple_cartesian<double> Kernel ;
typedef Polyhedron_3<Kernel,Polyhedron_items_with_id_3> Polyhedron ;

int sOK     = 0 ; 
int sFailed = 0 ;
 
#define CHECK(pred) \
        if (!(pred)) \
        { \
          cerr << "Assertion failure: " << #pred << endl \
               << "File:" << __FILE__ << endl \
               << "Line:" << __LINE__ << endl ; \
          throw 0 ; \
        }
        
#define CHECK_EQUAL(x,y)     CHECK(((x)==(y)))
#define CHECK_NOT_EQUAL(x,y) CHECK(((x)!=(y)))

