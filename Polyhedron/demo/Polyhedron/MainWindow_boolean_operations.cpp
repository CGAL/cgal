#include <QApplication>
#include <QTime>

#include "MainWindow.h"
#include "Scene.h"
#include "Polyhedron_type.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h> 
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Inverse_index.h>
#include <iostream>
#include <fstream>

// Boolean operations work only with exact kernel
typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_Kernel;
typedef CGAL::Polyhedron_3<Exact_Kernel> Exact_polyhedron;

// quick hacks to convert polyhedra from exact to inexact and vice-versa

template <class Polyhedron_input,
class Polyhedron_output>
struct Copy_polyhedron_to 
  : public CGAL::Modifier_base<typename Polyhedron_output::HalfedgeDS>
{
  Copy_polyhedron_to(const Polyhedron_input& in_poly) 
    : in_poly(in_poly) {}

  void operator()(typename Polyhedron_output::HalfedgeDS& out_hds)
  {
    typedef typename Polyhedron_output::HalfedgeDS Output_HDS;
    typedef typename Polyhedron_input::HalfedgeDS Input_HDS;

    CGAL::Polyhedron_incremental_builder_3<Output_HDS> builder(out_hds);

    typedef typename Polyhedron_input::Vertex_const_iterator Vertex_const_iterator;
    typedef typename Polyhedron_input::Facet_const_iterator  Facet_const_iterator;
    typedef typename Polyhedron_input::Halfedge_around_facet_const_circulator HFCC;

    builder.begin_surface(in_poly.size_of_vertices(),
      in_poly.size_of_facets(),
      in_poly.size_of_halfedges());


    for(Vertex_const_iterator
      vi = in_poly.vertices_begin(), end = in_poly.vertices_end();
      vi != end ; ++vi) 
    {
      typename Polyhedron_output::Point_3 p(::CGAL::to_double( vi->point().x()),
	::CGAL::to_double( vi->point().y()),
	::CGAL::to_double( vi->point().z()));
      builder.add_vertex(p);
    }

    typedef CGAL::Inverse_index<Vertex_const_iterator> Index;
    Index index( in_poly.vertices_begin(), in_poly.vertices_end());



    for(Facet_const_iterator 
      fi = in_poly.facets_begin(), end = in_poly.facets_end();
      fi != end; ++fi) 
    {
      HFCC hc = fi->facet_begin();
      HFCC hc_end = hc;
      //     std::size_t n = circulator_size( hc);
      //     CGAL_assertion( n >= 3);
      builder.begin_facet ();
      do {
	builder.add_vertex_to_facet(index[hc->vertex()]);
	++hc;
      } while( hc != hc_end);
      builder.end_facet();
    }
    builder.end_surface();
  } // end operator()(..)
private:
  const Polyhedron_input& in_poly;
}; // end Copy_polyhedron_to<>

template <class Poly_A, class Poly_B>
void copy_to(const Poly_A& poly_a, Poly_B& poly_b)
{
  Copy_polyhedron_to<Poly_A, Poly_B> modifier(poly_a);
  poly_b.delegate(modifier);
}

void from_exact(Exact_polyhedron& in,
		Polyhedron& out)
{
  copy_to(in, out);
  // 	std::ofstream out_stream("tmp.off");
  // 	out_stream << in;
  // 	std::ifstream in_stream("tmp.off");
  // 	in_stream >> out;
}

void to_exact(Polyhedron& in,
	      Exact_polyhedron& out)
{
  copy_to(in, out);
  // 	std::ofstream out_stream("tmp.off");
  // 	out_stream << in;
  // 	std::ifstream in_stream("tmp.off");
  // 	in_stream >> out;
  CGAL_assertion(out.is_valid());
}

void MainWindow::on_actionUnion_triggered()
{
  boolean_operation(BOOLEAN_UNION);
}

void MainWindow::on_actionIntersection_triggered()
{
  boolean_operation(BOOLEAN_INTERSECTION);
}

void MainWindow::on_actionDifference_triggered()
{
  boolean_operation(BOOLEAN_DIFFERENCE);
}

void MainWindow::boolean_operation(const Boolean_operation operation)
{
  const int indexA = scene->selectionAindex();
  const int indexB = scene->selectionBindex();

  Polyhedron* polyA = scene->polyhedron(indexA);
  Polyhedron* polyB = scene->polyhedron(indexB);
  if(!polyA) return;
  if(!polyB) return;
  if(polyA == polyB) return;

  QApplication::setOverrideCursor(Qt::WaitCursor);

  typedef CGAL::Nef_polyhedron_3<Exact_Kernel> Nef_polyhedron; 

  Exact_polyhedron exact_polyA; 
  to_exact(*polyA,exact_polyA);

  Exact_polyhedron exact_polyB; 
  to_exact(*polyB,exact_polyB);

  // convert to nef polyhedra
  QTime time;
  time.start();
  std::cout << "Convert to nef polyhedra...";
  Nef_polyhedron n1(exact_polyA); 
  Nef_polyhedron n2(exact_polyB); 
  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

  // perform Boolean operation
  std::cout << "Boolean operation...";
  time.start();
  switch(operation)
  {
  case BOOLEAN_UNION:
    n1 += n2;
    break;
  case BOOLEAN_INTERSECTION:
    n1 *= n2;
    break;
  case BOOLEAN_DIFFERENCE:
    n1 -= n2;
  }
  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

  // save the exact resulting mesh
  Exact_polyhedron exact_result;
  n1.convert_to_Polyhedron(exact_result);

  // reload as inexact one
  Polyhedron *pResult = new Polyhedron;
  from_exact(exact_result,*pResult);

  switch(operation)
  {
  case BOOLEAN_UNION:
    scene->addPolyhedron(pResult,
      tr("%1 union %2").arg(scene->polyhedronName(indexA),scene->polyhedronName(indexB)),
      Qt::yellow,
      true,
      scene->polyhedronRenderingMode(indexA));
    break;
  case BOOLEAN_INTERSECTION:
    scene->addPolyhedron(pResult,
      tr("%1 intersection %2").arg(scene->polyhedronName(indexA),scene->polyhedronName(indexB)),
      Qt::yellow,
      true,
      scene->polyhedronRenderingMode(indexA));
    break;
  case BOOLEAN_DIFFERENCE:
    scene->addPolyhedron(pResult,
      tr("%1 minus %2").arg(scene->polyhedronName(indexA),scene->polyhedronName(indexB)),
      Qt::yellow,
      true,
      scene->polyhedronRenderingMode(indexA));
  }

  QApplication::setOverrideCursor(Qt::ArrowCursor);
}

