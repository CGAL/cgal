// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/IO/PS_Stream_3.h
// package       : PS_Stream
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : 
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================
 
#ifndef CGAL_PS_STREAM_3_H
#define CGAL_PS_STREAM_3_H

#include <CGAL/Cartesian.h>
#include <cmath>
#include <list>
#include <CGAL/Bbox_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#include <CGAL/Polyhedron_iterator_3.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Triangulation_euclidean_traits_xy_3.h>

#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Regular_triangulation_2.h>

#include <CGAL/IO/Postscript_stream.h>

#include <CGAL/IO/PS_edge_3.h>
#include <CGAL/IO/PS_facet_3.h>
#include <CGAL/IO/PS_Stream.h>

CGAL_BEGIN_NAMESPACE

typedef double coord_type;
//typedef leda_real coord_type;

typedef CGAL::Cartesian< coord_type > CT;

typedef Bbox_3 PS_BBox3;

//typedef Cartesian <leda_real> R;
typedef Cartesian<double> D;

typedef Direction_3< D > Direction;
typedef Aff_transformation_3< CT > Transformation;
typedef Vector_3< CT > Vector3;
typedef Point_3< CT > Point3;
typedef Point_2< CT > Point2;
typedef Plane_3< CT > Plane3;

typedef CGAL::Polygon_traits_2< CT > Poly_Traits;
typedef Poly_Traits::Point_2 Point;
typedef std::list<Point> Container;
typedef CGAL::Polygon_2<Poly_Traits, Container> Polygon;
typedef CGAL::Triangle_3< D > Triangle3;
typedef CGAL::Segment_3< D > Segment3;
typedef CGAL::Tetrahedron_3< D > Tetrahedron;

// typedef Halfedge_data_structure_polyhedron_default_3<D> HDS;
// typedef Polyhedron_default_traits_3<D> Traits;
// typedef Polyhedron_3<Traits,HDS> Polyhedron;


// typedef Polyhedron::Facet Facet;
// typedef Polyhedron::Plane Plane;
// typedef Polyhedron::Halfedge_const_handle Halfedge_handle;
// typedef Polyhedron::Facet_const_iterator             FCI;
// typedef Polyhedron::Halfedge_around_facet_const_circulator HFCC;

//ADD FOR THE TRIANGULATION

//typedef CGAL::Triangulation_euclidean_traits_xy_3<Rp>  Gt;
//typedef CGAL::Triangulation_vertex_base_2<Gt> Vb; 
//typedef CGAL::Constrained_triangulation_face_base_2<Gt>   Fb; 
//typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
//typedef CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds> Triangulation;
//typedef Cartesian<double> Rp;

// typedef double Nt;
// typedef CGAL::Cartesian<Nt> Rp;
// typedef CGAL::Triangulation_euclidean_traits_xy_3<Rp> Gt;
// typedef Gt::Point    Point_tr ;
// typedef Gt::Segment  Segment_tr;
// typedef Gt::Triangle Triangle_tr;
// typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
// typedef CGAL::Constrained_triangulation_face_base_2<Gt>  Fb; //modif
// typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
// typedef CGAL::Constrained_Delaunay_triangulation_2<Gt, Tds> Triangulation;

// typedef Triangulation::Face_handle Face_handle;
// typedef Triangulation::Vertex_handle Vertex_handle;
// typedef Triangulation::Line_face_circulator  Line_face_circulator;
// typedef Triangulation::Face_iterator  Face_iterator;
//typedef Triangulation::Vertex_iterator  Vertex_iterator;


class PS_Stream_3 : public PS_Stream
{
public :

  // Constructors
  PS_Stream_3(const PS_BBox3& bb3,const Direction& d, const Direction &l,
              ostream& os,OutputMode = QUIET );
  PS_Stream_3(const PS_BBox3& bb3, const Direction& d, const Direction &l,
              const char* fname,OutputMode = QUIET);
  PS_Stream_3(const PS_BBox3& bb3,const Direction& d, const Direction &l,
              float H, ostream& os,OutputMode = QUIET );
  PS_Stream_3(const PS_BBox3& bb3,const Direction& d,const Direction &l,
              float H, const char* fname, OutputMode = QUIET);
  
  //Accessors 
  PS_BBox3    PS_Stream_3::bbox3() const {return _bbox3;}
  
  //Inline functions
  //Return the number of facet of the scene
  inline int number_of_face() const {return _number_of_face;}
  //Return the light
  inline Direction light() const { return _light;}
  //Return the point of view
  inline Direction point_of_view() const { return _point_of_view;}
  //Return the transformation of the scene to put the point_of_view
  //according the z-axis.
  Transformation transformation() const {return _transformation;}
  //Return the minimum in x of the array v
  coord_type search_xmin(Point3 v[]);
  //Return the maximum in x of the array v
  coord_type search_xmax(Point3 v[]); 
  //Return the minimum in y of the array v
  coord_type search_ymin(Point3 v[]); 
  //Return the maximum in y of the array v
  coord_type search_ymax(Point3 v[]);
  //Return the maximum in x,y,z
  coord_type search_max(PS_BBox3& b);
  //Return the minimum in x,y,z
  coord_type search_min(PS_BBox3& b);

  //Settings
  //Modify the type of filling of the facet
  void set_current_filling(const CGAL::FILLING& f){_current_filling=f;}
  //Modify the point of view vector
  void set_point_of_view(Direction& pov){_point_of_view=pov;}
  //Modify the light vector
  void set_light(Direction& light){_light=light;}

  //Operators
  //Return the norm
  coord_type norme(coord_type x,coord_type y,coord_type z){
    return sqrt(CGAL_NTS square(to_double(x))+
                CGAL_NTS square(to_double(y))+
                CGAL_NTS square(to_double(z)));}
  //Return the result of sqrt(1-x^2) used for the trtansformation matrix
  double den(double x){return (sqrt(1-x*x));}
  //Project a point3 and return the point2
  Point2 transform( const Point3& p);
  //Transform a point3 according the transformation and return the point2
  Point2 transform( Transformation t, const Point3& p); 
  //This method allows to convert a point in Leda_Real in a point2 in double
  Point_2<D> PS_Stream_3::convert_to_double(Point2& p);
  ///**************************************************************//
  //Make the sort according the increasing z
  vector<PS_facet_3> PS_Stream_3::zmin_sort();
  //make the normal sort to determine visible surfaces
  vector<PS_facet_3> PS_Stream_3::normal_sort();
  //Display the scene
  void display();
  //Transform all the facets before the display of the scene.
  // These transformation allows to
  //transform the point of view to become along the z-axis.
  void PS_Stream_3::transformation();
  //Allows to split the face number i according the intersection between the
  //face i and j
  int PS_Stream_3::cutting(vector<PS_facet_3>& liste_face,int i,int j);
  //Make the depth sort
  void PS_Stream_3::depth_sort();
  //Check if the face i is under the face j
  int PS_Stream_3::dessous(PS_facet_3& facei,PS_facet_3& facej );


  //Adding functions 
  //Allows to add a new vector of PS_facet_3 in the PS_Stream_3
  void PS_Stream_3::add_vfacet(vector<PS_facet_3> &vfacet);
  //Allows to add a new  PS_facet_3 in the PS_Stream_3
  void PS_Stream_3::add_facet(PS_facet_3 &face);
  //Allows to add a point 3D in the PS_Stream_3
  void PS_Stream_3::add_point(Point3 &point);
  //Allows to add a segment 3D in the PS_Stream_3
  void PS_Stream_3::add_segment(Segment3 &segment);
  //Allows to add a Line 3D in the PS_Stream_3
  void PS_Stream_3::add_line(Line3 &line);
  //Allows to add a Triangulation_2 in the PS_Stream_3
  template <class triangulation_2>
  void PS_Stream_3::add_triangulation(triangulation_2& tr);
  //Allows to add a triangle 3D in the PS_Stream_3
  void PS_Stream_3::add_triangle(Triangle3 &triangle);
  //Allows to add a tetrahedron in the PS_Stream_3
  void PS_Stream_3::add_tetrahedron(Tetrahedron &t);
  //Allows to add a polyhedron in the PS_Stream_3
  template <class Traits, class HDS>
  void PS_Stream_3::add_polyhedron(Polyhedron_3<Traits,HDS>& p);

  //This method transform a point 3D in a vector of PS_facet_3..}
  //It is necessary to apply the algorithm to a point.
  //So the point is assimilated to a little cube
  vector<PS_facet_3> point_to_facet(Point3& point);
   //This method transform a segment 3D in a vector of PS_facet_3.
  //It is necessary to apply the algorithm to a point.
  //So the segment is assimilated to a parrallepipede
  vector<PS_facet_3> segment_to_facet(Segment3& segment);
  //This method transform a line 3D in a vector of PS_facet_3.
  //It is necessary to apply the algorithm to a point.
  //So the segment is assimilated
  //to a parrallepipede which is bigger than the BBox
  vector<PS_facet_3> line_to_facet(Line3& line);
  //This method transform a tetrahedron into a vector of PS_facet_3.
  vector<PS_facet_3> tetrahedron_to_facet(Tetrahedron& T1);
  //This method transform a triangle 3D in PS_facet_3.
  PS_facet_3 triangle_to_facet(Triangle3& T);
  //This method transform a triangulation 2 into a vector of PS_facet_3.
  template <class triangulation_2>
  vector<PS_facet_3> triangulation_2_to_facet(triangulation_2& tr);
  //This method transform a polyhedron_3 into a vector of PS_facet_3.
  template <class Traits, class HDS>
  vector<PS_facet_3> polyhedron_to_facet(Polyhedron_3<Traits,HDS>& p); 
  
  //Output
  //Allows to add a tetrahedron in the PS_Stream_3
  friend PS_Stream_3& operator << (PS_Stream_3& ps, Tetrahedron& t);
  //Allows to add a point 3D in the PS_Stream_3
  friend PS_Stream_3& operator << (PS_Stream_3& ps,  Point3& p);
  //Allows to add a PS_edge_3 in the PS_Stream_3
  friend PS_Stream_3& operator << (PS_Stream_3& ps, const PS_edge_3& ar);
  //Allows to add a PS_facet_3 in the PS_Stream_3
  friend PS_Stream_3& operator << (PS_Stream_3& ps, const PS_facet_3& face);
  //Allows to add a Segment_3 in the PS_Stream_3
  friend PS_Stream_3& operator << (PS_Stream_3& ps, Segment3& segment);
  //Allows to add a Line_3 in the PS_Stream_3
  friend PS_Stream_3& operator << (PS_Stream_3& ps, Line3& line);
  //Allows to add a Triangle_3 in the PS_Stream_3
  friend PS_Stream_3& operator << (PS_Stream_3& ps,  Triangle3& t);
  //Allows to modify the filling type of the PS_Stream_3
  friend PS_Stream_3& operator << (PS_Stream_3& ps, CGAL::FILLING& filling);
  //Allows to add a Polyhedron_3 in the PS_Stream_3
  template <class Traits, class HDS>
  friend PS_Stream_3& operator <<(PS_Stream_3& ps,Polyhedron_3<Traits,HDS>& p);
  //Allows to add a Triangulation_2 in the PS_Stream_3
  template <class Gt,class Tds>
  friend PS_Stream_3& operator << (PS_Stream_3& ps,Triangulation_2<Gt,Tds> &t);
  //Allows to add a Delaunay_triangulation_2 in the PS_Stream_3
  template < class Gt, class Tds >
  friend PS_Stream_3& operator << (PS_Stream_3& ps,
                                   Delaunay_triangulation_2<Gt,Tds> &t);
  //Allows to add a Constrained_triangulation_2 in the PS_Stream_3
  template < class Gt, class Tds>
  friend PS_Stream_3& operator<<(PS_Stream_3& ps,
                                 Constrained_triangulation_2<Gt,Tds> &t);
  //Allows to add a Regular_triangulation_2 in the PS_Stream_3
  template < class Gt, class Tds >
  friend PS_Stream_3& operator << (PS_Stream_3& ps,
                                   Regular_triangulation_2<Gt,Tds> &t); 
  //Allows to add a Constrained_Delaunay_triangulation_2 in the PS_Stream_3
  template < class Gt, class Tds >
  friend PS_Stream_3& operator << (PS_Stream_3& ps,
                              Constrained_Delaunay_triangulation_2<Gt,Tds> &t);
  //AJOUT MANIPULATOR
  //Modify the type of filling of the facet
  //friend PS_Stream_3& set_point_of_view(Direction& pov);
    //{_point_of_view=pov;}
  //Modify the color of the facet
  //friend PS_Stream_3& set_light(Direction& light);// {_light=light;}

  //  PS_Manipulator_creator<const
  //extern PS_Manipulator_creator<const Direction&>
  //p_o_v(&PS_Stream_3::set_point_of_view);
  //extern PS_Manipulator_creator<const Direction&>
  //_light_(&PS_Stream_3::set_light);
  //The BBox3
  // Define the bounding box
  PS_BBox3 _bbox3;

private :

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
  
  //Define the vector of PS_facet_3 of the PS_Stream_3
  vector<PS_facet_3> _scene_initiale;
  //Define the vector of visible PS_facet_3 of the PS_Stream_3
  //This allows for the user to not re-enter all the facet if the
  //scene point of view has changed. So The display will only work
  //with this vector.
  vector<PS_facet_3> _scene;
   //number of facet
  int _number_of_face;
  // Define the point of view direction
  Direction _point_of_view;
  // Define the light direction
  Direction _light;
  //Define the transformation matrice
  Transformation _transformation;
  //Context
  //The current style of filling
  FILLING _current_filling;
};

CGAL_END_NAMESPACE

#endif
