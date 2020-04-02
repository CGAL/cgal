#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

template <class Kernel, class HDS>
class gausian_map_to_polyhedron_3 : public CGAL::Modifier_base<HDS> {

  typedef CGAL::Gausian_map<Kernel> Gausian_map;
  typedef typename Gausian_map::SFace_const_iterator SFace_const_iterator;
  typedef typename Gausian_map::SVertex_const_iterator SVertex_const_iterator;
  typedef typename Gausian_map::SHalfedge_around_svertex_const_circulator
    SHalfedge_around_svertex_const_circulator;

  const Gausian_map& G;
  CGAL::Unique_hash_map<SFace_const_iterator, int> SFace2int;

 public:
  gausian_map_to_polyhedron_3(const Gausian_map& Gin) : G(Gin) {}

    void operator()( HDS& hds) {
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( G.number_of_sfaces(), G.number_of_svertices(), G.number_of_shalfedges());
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;

        int i = 0;
        SFace_const_iterator fi;
        for(fi = G.sfaces_begin(); fi != G.sfaces_end(); ++fi) {
          B.add_vertex(fi->mark().point());
          SFace2int[fi] = i++;
        }

        SVertex_const_iterator vi;
        for(vi = G.svertices_begin(); vi != G.svertices_end(); ++vi) {
          SHalfedge_around_svertex_const_circulator
            svc(G.first_out_edge(vi)),
            svend(svc);
          B.begin_facet();
          CGAL_For_all(svc,svend)
            B.add_vertex_to_facet(SFace2int[svc->incident_sface()]);
          B.end_facet();
        }
        B.end_surface();
    }
};
