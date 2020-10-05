#ifndef CGAL_NEF3_SAT_WRITER
#define CGAL_NEF3_SAT_WRITER

#include <CGAL/Nef_3/SNC_decorator.h>

namespace CGAL {

template <typename Nef_>
class sat_writer : public CGAL::SNC_decorator<typename Nef_::SNC_structure> {
  typedef Nef_ Nef_polyhedron;
  typedef typename Nef_::SNC_structure        SNC_structure;
  typedef CGAL::SNC_decorator<SNC_structure>  Base;

  typedef typename SNC_structure::Volume_iterator Volume_iterator;
  typedef typename SNC_structure::Vertex_iterator Vertex_iterator;
  typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator;
  typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator;
  typedef typename SNC_structure::SHalfedge_iterator SHalfedge_iterator;
  typedef typename SNC_structure::SHalfloop_iterator SHalfloop_iterator;
  typedef typename SNC_structure::SFace_iterator SFace_iterator;

  typedef typename SNC_structure::Volume_handle Volume_handle;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle SFace_handle;

  typedef typename SNC_structure::Plane_3 Plane_3;
  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Vector_3 Vector_3;

 private:
  std::ostream& out;

  CGAL::Unique_hash_map<Vertex_iterator,int> VI;
  CGAL::Unique_hash_map<Halfedge_iterator,int> EI;
  CGAL::Unique_hash_map<Halffacet_iterator,int>    FI;
  CGAL::Unique_hash_map<Volume_iterator,int>   CI;
  CGAL::Unique_hash_map<SHalfedge_iterator,int> SEI;
  CGAL::Unique_hash_map<SHalfloop_iterator,int>   SLI;
  CGAL::Unique_hash_map<SFace_iterator,int>     SFI;

  CGAL::Unique_hash_map<Vertex_iterator,int> VERTEX;
  CGAL::Unique_hash_map<Halffacet_iterator,int> FACE;
  CGAL::Unique_hash_map<Volume_iterator,int> LUMP;
  CGAL::Unique_hash_map<Halfedge_iterator,int> COEDGE;
  CGAL::Unique_hash_map<Halfedge_iterator,int> EDGE;
  CGAL::Unique_hash_map<Halfedge_iterator,int> NEXT_COEDGE;
  CGAL::Unique_hash_map<Halfedge_iterator,int> PREV_COEDGE;

  CGAL::Unique_hash_map<Halffacet_iterator, int> NEXT_FACE;
  CGAL::Unique_hash_map<Volume_iterator,int> FIRST_FACE;

  int point_offset;
  int edge_offset;
  int coedge_offset;
  int face_offset;
  int loop_offset;

 public:
  sat_writer(std::ostream& os, Nef_polyhedron& N) :
    Base(*const_cast<SNC_structure*>(N.sncp())), out(os) {

    Vertex_iterator v=vertices_begin();
    for(int i=0;v!=vertices_end();++v)
      VI[v]=i++;
    Halfedge_iterator e=halfedges_begin();
    for(int i=0;e!=halfedges_end();++e)
      EI[e]=i++;
    Halffacet_iterator f=halffacets_begin();
    for(int i=0;f!=halffacets_end();++f)
      FI[f]=i++;
    Volume_iterator c=volumes_begin();
    for(int i=0;c!=volumes_end();++c)
      CI[c]=i++;
    SHalfedge_iterator se=shalfedges_begin();
    for(int i=0;se!=shalfedges_end();++se)
      SEI[se]=i++;
    SHalfloop_iterator sl=shalfloops_begin();
    for(int i=0;sl!=shalfloops_end();++sl)
      SLI[sl]=i++;
    SFace_iterator sf=sfaces_begin();
    for(int i=0;sf!=sfaces_end();++sf)
      SFI[sf]=i++;

    int i=2;
    for(c=++volumes_begin();c!=volumes_end();++c) {
      LUMP[c]=i;
      i+=2;
    }

    for(f=halffacets_begin();f!=halffacets_end();++f)
      if(CI[f->incident_volume()]==0) {
        FACE[f]=i;
        i+=3;
      }

    for(se=shalfedges_begin();se!=shalfedges_end();++se)
      if(CI[se->facet()->incident_volume()]==0)
        COEDGE[se->source()]=i++;

    for(se=shalfedges_begin();se!=shalfedges_end();++se)
      if(CI[se->facet()->incident_volume()]==0) {
        NEXT_COEDGE[se->source()] = COEDGE[se->next()->source()];
        PREV_COEDGE[se->source()] = COEDGE[se->prev()->source()];
      }

    for(se=shalfedges_begin();se!=shalfedges_end();++se)
      if(CI[se->facet()->incident_volume()]==0 && se->source()->is_twin()) {
        EDGE[se->source()]=i;
        i+=2;
      }

    for(v=vertices_begin();v!=vertices_end();++v) {
      VERTEX[v]=i;
      i+=2;
    }

    std::list<Halffacet_handle> next_face[number_of_volumes()-1];
    for(f=halffacets_begin();f!=halffacets_end();++f)
      if(CI[f->incident_volume()]==0)
        next_face[CI[f->twin()->incident_volume()]-1].push_back(f);

    for(c=++volumes_begin();c!=volumes_end();++c)
      FIRST_FACE[c] = FACE[*next_face[CI[c]-1].begin()];

    typename std::list<Halffacet_handle>::const_iterator li;
    for(i=0; i<number_of_volumes()-1;++i) {
      for(li=next_face[i].begin();li!=--next_face[i].end();)
        NEXT_FACE[*li] = FACE[*++li];
      NEXT_FACE[*li] = -1;
    }
  }

  void print_header() {
    out << "700 0 1 0" << std::endl
        << "22 ACIS/Scheme AIDE - 7.0 11 ACIS 7.0 NT 24 Mon Oct 18 18:53:23 2004" << std::endl
        << "1 9.9999999999999995e-007 1e-010" << std::endl;
  }

  void print_footer() {
    out << "End-of-ACIS-data" << std::endl;
  }

  void print_body() {
    out << "body $-1 -1 $-1 $2 $-1 $1 #" << std::endl;
    out << "transform $-1 -1 1 0 0 0 1 0 0 0 1 0 0 0 1 no_rotate no_reflect no_shear #"
        << std::endl;
  }

  void print_lumps_and_shells() {

    int i=0;
    for(Volume_iterator c=++volumes_begin(); c!=volumes_end(); ++c) {

      Volume_iterator c_next(c);
      ++c_next;
      SFace_handle sf(c->shells_begin());
      SHalfedge_handle e(sf->sface_cycles_begin());
      Halffacet_handle f(e->twin()->facet());
      if(CI[f->incident_volume()]!=0)
        f=f->twin();

      int next = (c_next!=volumes_end()?LUMP[c_next]:-1);
      out << "lump $-1 -1 $-1 $" << next
          << " $" << 3+2*i << " $0 #" << std::endl;
      out << "shell $-1 -1 $-1 $-1" // << (next==-1?-1:next+1)
          << " $-1 $" << FIRST_FACE[c]
          << " $-1 $" << 2+2*i << " #" << std::endl;
      ++i;
    }
  }

  void print_face() {

    Halffacet_iterator f;
    for(f=halffacets_begin();f!=halffacets_end();++f) {
      if(CI[f->incident_volume()]!=0) continue;

      out << "face $-1 -1 $-1 $" << NEXT_FACE[f] << " $" << FACE[f]+2 << " $"
          << LUMP[f->twin()->incident_volume()]+1 << " $-1 $"
          << FACE[f]+1 << " forward single #" << std::endl;

      Plane_3 pl(f->plane());
      SHalfedge_handle se(f->facet_cycles_begin());
      Halfedge_handle e(se->source());
      Point_3 pt(e->source()->point());
      Vector_3 vec(e->point()-CGAL::ORIGIN);

      out << "plane-surface $-1 -1 $-1 "
          << CGAL::to_double(pt.x()) << " "
          << CGAL::to_double(pt.y()) << " "
          << CGAL::to_double(pt.z()) << " "
          << pl.a() << " "
          << pl.b() << " "
          << pl.c() << " "
          << vec.hx() << " "
          << vec.hy() << " "
          << vec.hz() << " "
          << "forward_v I I I I #" << std::endl;

      out << "loop $-1 -1 $-1 $-1 $"
          << COEDGE[se->source()] << " $" << FACE[f] << " #" << std::endl;
    }
  }

  void print_coedge() {

    SHalfedge_iterator se;
    for(se=shalfedges_begin();se!=shalfedges_end();++se)
      if(CI[se->facet()->incident_volume()]==0) {
        out << "coedge $-1 -1 $-1 $" << NEXT_COEDGE[se->source()]
            << " $" << PREV_COEDGE[se->source()]
            << " $" << COEDGE[se->source()->twin()] << " $";
        if(se->source()->is_twin())
          out << EDGE[se->source()] << " reversed $";
        else
          out << EDGE[se->source()->twin()] << " forward $";
        out << FACE[se->facet()]+2
            << " $-1 #" << std::endl;
      }
  }

  void print_edges() {
    SHalfedge_iterator se;
    for(se=shalfedges_begin();se!=shalfedges_end();++se)
      if(CI[se->facet()->incident_volume()]==0 && se->source()->is_twin()) {

        double end_param;
        Vector_3 vec(se->source()->twin()->source()->point() -
                     se->source()->source()->point());
        if(vec.x()!=0)
          end_param = CGAL::to_double(vec.x()/se->source()->point().x());
        else if(vec.y()!=0)
          end_param = CGAL::to_double(vec.y()/se->source()->point().y());
        else
          end_param = CGAL::to_double(vec.z()/se->source()->point().z());

        out << "edge $-1 -1 $-1 $"
            << VERTEX[se->source()->source()] << " 0 $"
            << VERTEX[se->source()->twin()->source()] << " " << end_param << " $"
            << COEDGE[se->source()]
            << " $" << EDGE[se->source()]+1
            << " forward @7 unknown #" << std::endl;

        out << "straight-curve $-1 -1 $-1 "
            << CGAL::to_double(se->source()->source()->point().x()) << " "
            << CGAL::to_double(se->source()->source()->point().y()) << " "
            << CGAL::to_double(se->source()->source()->point().z()) << " "
            << CGAL::to_double(se->source()->point().x()) << " "
            << CGAL::to_double(se->source()->point().y()) << " "
            << CGAL::to_double(se->source()->point().z()) << " I I #"
            << std::endl;
      }
  }

  void print_vertex() {
    Vertex_iterator v;
    for(v=vertices_begin();v!=vertices_end();++v) {
      Halfedge_handle e=v->svertices_begin();
      if(!e->is_twin()) e=e->twin();
      out << "vertex $-1 -1 $-1 $" << EDGE[e]
          << " $" << VERTEX[v]+1 << " #" << std::endl;

      out << "point $-1 -1 $-1 "
          << CGAL::to_double(v->point().x()) << " "
          << CGAL::to_double(v->point().y()) << " "
          << CGAL::to_double(v->point().z()) << " #"
          << std::endl;

    }
  }

  void print() {
    print_header();
    print_body();
    print_lumps_and_shells();
    print_face();
    print_coedge();
    print_edges();
    print_vertex();
    print_footer();
  }
};

} //namespace CGAL
#endif // CGAL_NEF3_SAT_WRITER
