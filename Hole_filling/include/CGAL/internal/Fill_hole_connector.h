#ifndef CGAL_HOLE_FILLING_FILL_HOLE_CONNECTOR
#define CGAL_HOLE_FILLING_FILL_HOLE_CONNECTOR
#include <CGAL/internal/Fair_Polyhedron_3.h>
#include <CGAL/internal/Refine_Polyhedron_3.h>
#include <CGAL/internal/Triangulate_hole_Polyhedron_3.h>
#include <CGAL/trace.h>

#include <boost/function_output_iterator.hpp>

#include <vector>
#include <set>
#include <map>

namespace CGAL {
namespace internal {

  struct Nop_functor {
    template<class T>
    void operator()(const T& t) const {}
  };

  // Connector class to triangulate - refine - fair holes on Polyhedron_3
  template<class Polyhedron>
  class Fill_hole_Polyhedron_3 {

    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Facet_handle Facet_handle;
    typedef typename Polyhedron::Halfedge_around_facet_circulator  Halfedge_around_facet_circulator;
    typedef typename Polyhedron::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;

  private:
    void average_length(Polyhedron& poly, Vertex_handle vh, std::map<Vertex_handle, double>& scale_attribute)
    {
      const Point_3& vp = vh->point(); 
      Halfedge_around_vertex_circulator circ(vh->vertex_begin()), done(circ);
      int deg = 0;
      double sum = 0;
      do {
        const Point_3& vq = circ->opposite()->vertex()->point();
        sum += std::sqrt(CGAL::squared_distance(vp, vq));
        ++deg;
        ++circ;
      } while(circ != done);
      scale_attribute[vh] = sum/deg;
    }

    void get_boundary_vertices(Polyhedron& poly, Halfedge_handle it, std::set<Vertex_handle>& boundary_vertices) {
      Halfedge_around_facet_circulator circ(it), done(circ);
      do{
        boundary_vertices.insert(circ->vertex()); 
      } while (++circ != done);
    }

  public:
    template<class OutputIterator>
    void triangulate_and_refine_hole(Polyhedron& poly, Halfedge_handle it, double alpha, OutputIterator output)
    {
      if(! it->is_border()){ return; }

      // compute scale_attribute before triangulation
      std::map<Vertex_handle, double> scale_attribute;
      Halfedge_around_facet_circulator circ(it), done(circ);
      do{
        average_length(poly, circ->vertex(), scale_attribute); 
      } while (++circ != done);
      // Triangulate
      //std::set<Facet_handle> facets;
      
      //std::cout << "is_valid before: " << poly.is_valid() << std::endl;
      std::vector<Facet_handle> facets;
      triangulate_hole(poly, it, std::back_inserter(facets));//std::inserter(facets, facets.begin()));
      //std::cout << "is_valid after: " << poly.is_valid() << std::endl;
      //std::cout << "f " << facets.size() << std::endl;
      // Refine
      internal::Refine_Polyhedron_3<Polyhedron> refine_functor(alpha);
      refine_functor(poly, scale_attribute, facets.begin(), facets.end(), output);
    }

    template<class SparseLinearSolver, class WeightCalculator>
    void triangulate_refine_and_fair
      (Polyhedron& poly, Halfedge_handle it, double alpha, WeightCalculator weight_calculator) 
    {
      if(! it->is_border()){ return; }

      // save boundary vertices before triangulation
      std::set<Vertex_handle> boundary_vertices;
      get_boundary_vertices(poly, it, boundary_vertices);

      //Triangulate and refine
      std::vector<Facet_handle> facets;
      triangulate_and_refine_hole(poly, it, alpha, std::back_inserter(facets));

      // get interior vertices
      std::set<Vertex_handle> interior_vertices;
      for(typename std::vector<Facet_handle>::iterator it = facets.begin(); it != facets.end(); ++it) {
        Halfedge_around_facet_circulator circ = (*it)->facet_begin();
        do {
          if(boundary_vertices.find(circ->vertex()) == boundary_vertices.end()) {
            interior_vertices.insert(circ->vertex());
          }
        } while(++circ != (*it)->facet_begin());
      }

      CGAL_TRACE_STREAM << "before fair |boundary vertices| = " << boundary_vertices.size() << std::endl;
      CGAL_TRACE_STREAM << "before fair |interior vertices| = " << interior_vertices.size() << std::endl;
      // Fair
      fair<SparseLinearSolver>(poly, interior_vertices, weight_calculator);
    }
  };
}//namespace internal

/**
 * @brief Function triangulating a hole in surface mesh.
 * 
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam OutputIterator iterator holding 'Polyhedron::Facet_handle' for patch facets.
 * @param[in, out] polyhedron surface mesh which has the hole
 * @param border_halfedge a border halfedge incident to the hole
 * @param[out] output iterator over patch facets. It can be omitted if output is not required.
 * 
 * \warning Using this function on very large holes might not be feasible, since the cost of triangulation is O(n^3).
 */
template<class Polyhedron, class OutputIterator>
void triangulate_hole(Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  OutputIterator output
  ) 
{
  internal::Triangulate_hole_Polyhedron_3<Polyhedron>()(polyhedron, border_halfedge, output);
}

template<class Polyhedron>
void triangulate_hole(Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge
  ) 
{
  triangulate_hole(polyhedron, border_halfedge, 
    boost::make_function_output_iterator(CGAL::internal::Nop_functor()));
}

}//namespace CGAL
#endif // CGAL_HOLE_FILLING_FILL_HOLE_CONNECTOR