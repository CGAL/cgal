#pragma warning(disable: 4267)

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h> 
#include <CGAL/Simple_cartesian.h>
#include <CGAL/tags.h> 


#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Inventor_ostream.h>
#include <CGAL/IO/VRML_1_ostream.h>
#include <CGAL/IO/VRML_2_ostream.h>
#include <CGAL/IO/Geomview_stream.h>

#include <CGAL/IO/File_writer_OFF.h>
#include <CGAL/IO/File_writer_inventor.h>
#include <CGAL/IO/File_writer_VRML_2.h>



#include <list>
#include <iostream>
#include <istream>
#include <fstream>


typedef CGAL::Simple_cartesian<double>      Kernel;
typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;
typedef Polyhedron::Facet                   Facet;
typedef Polyhedron::Facet_handle            Facet_handle;
typedef Polyhedron::Facet_const_handle      Facet_const_handle;
typedef Polyhedron::Halfedge_around_facet_circulator
                                            Halfedge_around_facet_circulator;
typedef Polyhedron::Halfedge_handle         Halfedge_handle;
typedef Polyhedron::Halfedge_const_handle   Halfedge_const_handle;
//typedef Polyhedron::Supports_removal        Supports_removal;
//typedef Polyhedron::Supports_vertex_halfedge Supports_vertex_halfedge;

typedef CGAL::Plane_3<Kernel>               Plane_3;
typedef CGAL::Point_3<Kernel>               Point_3;
typedef CGAL::Vector_3<Kernel>              Vector_3;
typedef CGAL::Tetrahedron_3<Kernel>         Tetrahedron_3;
typedef CGAL::Aff_transformation_3<Kernel>  Aff;

