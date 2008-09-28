#include <QApplication>
#include <QTime>

#include "MainWindow.h"
#include "Scene.h"
#include "Polyhedron_type.h"
#include "Nef_type.h"

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Inverse_index.h>
#include <iostream>
#include <fstream>

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
  CGAL_assertion(out.is_valid());
}

void to_exact(Polyhedron& in,
	      Exact_polyhedron& out)
{
  copy_to(in, out);
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

  if(indexA == indexB) return;

  Polyhedron* polyA = scene->polyhedron(indexA);
  Polyhedron* polyB = scene->polyhedron(indexB);
  Nef_polyhedron* nefpolyA = scene->nefPolyhedron(indexA);
  Nef_polyhedron* nefpolyB = scene->nefPolyhedron(indexB);

  if(!polyA && !nefpolyA) return;
  if(!polyB && !nefpolyB) return;

  QApplication::setOverrideCursor(Qt::WaitCursor);

  Nef_polyhedron* n1;
  Nef_polyhedron* n2;

  if(nefpolyA) {
    n1 = new Nef_polyhedron(*nefpolyA);
  }
  else {
    QTime time;
    time.start();
    std::cout << "Convert polyhedron A to nef polyhedron...";
    Exact_polyhedron exact_polyA; 
    to_exact(*polyA,exact_polyA);
    n1 = new Nef_polyhedron(exact_polyA); 
    std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
  }

  if(nefpolyB) {
    n2 = nefpolyB; // no copy, because n2 is not modified
  }
  else {
    QTime time;
    time.start();
    std::cout << "Convert polyhedron B to nef polyhedron...";
    Exact_polyhedron exact_polyB; 
    to_exact(*polyB, exact_polyB);
    n2 = new Nef_polyhedron(exact_polyB); 
    std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
  }
  // perform Boolean operation
  std::cout << "Boolean operation...";
  QTime time;
  time.start();
  switch(operation)
  {
  case BOOLEAN_UNION:
    (*n1) += (*n2);
    break;
  case BOOLEAN_INTERSECTION:
    (*n1) *= (*n2);
    break;
  case BOOLEAN_DIFFERENCE:
    (*n1) -= (*n2);
  }
  std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

  QString name;
  switch(operation)
  {
  case BOOLEAN_UNION:
    name = tr("%1 union %2");
    break;
  case BOOLEAN_INTERSECTION:
    name = tr("%1 intersection %2");
    break;
  case BOOLEAN_DIFFERENCE:
    name = tr("%1 minus %2");
  }
  
  if(n1->is_simple())
  {
    std::cout << "Convert result back to polyhedron...";
    time.restart();
    // save the exact resulting mesh
    Exact_polyhedron exact_result;
    n1->convert_to_Polyhedron(exact_result);
    // reload as inexact one
    Polyhedron *pResult = new Polyhedron;
    from_exact(exact_result,*pResult);
    std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;
    scene->addPolyhedron(pResult,
      name.arg(scene->polyhedronName(indexA),scene->polyhedronName(indexB)),
      Qt::magenta,
      true,
      scene->polyhedronRenderingMode(indexA));
  }
  else {
    scene->addNefPolyhedron(n1,
      name.arg(scene->polyhedronName(indexA),scene->polyhedronName(indexB)),
      Qt::green,
      true,
      scene->polyhedronRenderingMode(indexA));
  }

  if(!nefpolyB)
    delete n2; // destroyed n2, if it was a copy

  QApplication::restoreOverrideCursor();
}

