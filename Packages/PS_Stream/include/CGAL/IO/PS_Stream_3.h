#ifndef Postscript_STREAM_3
#define Postscript_STREAM_3

#include <math.h>
#include <list>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/Direction_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#include <CGAL/Polyhedron_iterator_3.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Polygon_2.h>

#include <CGAL/IO/Postscript_stream.h>

CGAL_BEGIN_NAMESPACE

typedef Bbox_3 PS_BBox3;
typedef Direction_3< Cartesian<double> > Direction;
typedef Aff_transformation_3< Cartesian<double> > Transformation;
typedef Vector_3< Cartesian<double> > Vector3;
typedef Point_3< Cartesian<double> > Point3;
typedef Point_2< Cartesian<double> > Point2;
typedef Cartesian<double> R;

typedef CGAL::Polygon_traits_2<R> Poly_Traits;
typedef Poly_Traits::Point_2 Point;
typedef std::list<Point> Container;
typedef CGAL::Polygon_2<Poly_Traits, Container> Polygon;



typedef Halfedge_data_structure_polyhedron_default_3<R> HDS;
typedef Polyhedron_default_traits_3<R> Traits;
typedef Polyhedron_3<Traits,HDS> Polyhedron;
typedef Polyhedron::Facet Facet;
typedef Polyhedron::Plane Plane;
typedef Polyhedron::Halfedge_const_handle Halfedge_handle;
typedef Polyhedron::Facet_const_iterator             FCI;
typedef Polyhedron::Halfedge_around_facet_const_circulator HFCC;


class PS_Stream_3 : public PS_Stream
{

private :

  // Define the bounding box
PS_BBox3 _bbox3;

  // Define the point of view direction
Direction _dir;

  //Define the transformation matrices
Transformation _t;

  //Array of vertices of the bounding box
 Point3 v[8];

  //Initialisation du tableau contenant les sommets de la bounding box
  void init(){ v[0]=Point3(_bbox3.xmin(),_bbox3.ymin(),_bbox3.zmin());
               v[1]=Point3(_bbox3.xmin(),_bbox3.ymin(),_bbox3.zmax());
               v[2]=Point3(_bbox3.xmin(),_bbox3.ymax(),_bbox3.zmin());
               v[3]=Point3(_bbox3.xmin(),_bbox3.ymax(),_bbox3.zmax());
               v[4]=Point3(_bbox3.xmax(),_bbox3.ymax(),_bbox3.zmax());
               v[5]=Point3(_bbox3.xmax(),_bbox3.ymin(),_bbox3.zmin());
               v[6]=Point3(_bbox3.xmax(),_bbox3.ymin(),_bbox3.zmax());
               v[7]=Point3(_bbox3.xmax(),_bbox3.ymax(),_bbox3.zmin());
             }



public :


  // Constructor

PS_Stream_3(const PS_BBox3& bb3,const Direction& d, ostream& os,
	    OutputMode = QUIET );
PS_Stream_3(const PS_BBox3& bb3, const Direction& d, const char* fname, 
	    OutputMode = QUIET);
PS_Stream_3(const PS_BBox3& bb3, const Direction& d, float H, const
	    char* fname, OutputMode = QUIET);
PS_Stream_3(const PS_BBox3& bb3,const Direction& d, float H, ostream& os,
	    OutputMode = QUIET );

  //Accessors
  PS_BBox3    bbox3() const {return _bbox3;}
  Direction  direction() const {return _dir;}
  Transformation trans() const {return _t;}


  //Utils
  double norme(double x,double y,double z){return
					     (sqrt(x*x+y*y+z*z));}
  double den(double x){return (sqrt(1-x*x));}


  Point2 transform( Transformation t, const Point3& p);
  double search_xmin(Point3 v[]);
  double search_xmax(Point3 v[]);
  double search_ymin(Point3 v[]);
  double search_ymax(Point3 v[]);
  Plane compute_plane_equations(const Facet& f);

};

template < class R >
PS_Stream_3& operator <<(PS_Stream_3& ps, Point_3<R>& p)
{
Point2 p2;
p2=ps.transform(ps.trans(),p);
ps << p2;
return ps;
}



CGAL_END_NAMESPACE

#endif // Postscript_STREAM_3


  
