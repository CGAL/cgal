#ifndef CGAL_VORONOI_DIAGRAM_2_FINDER_CLASSES_H
#define CGAL_VORONOI_DIAGRAM_2_FINDER_CLASSES_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>
#include <map>
#include <CGAL/Triangulation_utils_2.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//-------------------------------------------------------------------
//-------------------------------------------------------------------


template<class VDA>
struct Find_next_halfedge
{
  typedef Triangulation_cw_ccw_2  CW_CCW_2;

  typedef typename VDA::Dual_face_handle         Dual_face_handle;
  typedef typename VDA::Edge_degeneracy_tester   Edge_tester;

  void operator()(const VDA* vda, const Dual_face_handle& f, int i,
		  Dual_face_handle& fnext, int& inext) const
  {
    Dual_face_handle fcur = f;
    int icur = i, cw_i;
    do {
      cw_i = CW_CCW_2::cw(icur);
      fnext = fcur->neighbor(cw_i);
      inext = vda->dual_graph_data_structure().mirror_index(fcur, cw_i);
      fcur = fnext;
      icur = inext;
    } while ( vda->edge_tester()(fnext, inext) );
  }
};

//-------------------------------------------------------------------

template<class VDA>
struct Find_opposite_halfedge
{
  typedef typename VDA::Dual_face_handle     Dual_face_handle;
  typedef Find_next_halfedge<VDA>            Find_next_halfedge;

  void operator()(const VDA* vda, const Dual_face_handle& f, int i,
		  Dual_face_handle& fopp, int& iopp) const 
  {
    Dual_face_handle f1;
    int i1;

    int i_mirror = vda->dual_graph_data_structure().mirror_index(f, i);
    Find_next_halfedge()(vda, f->neighbor(i), i_mirror, f1, i1);

    fopp = f1->neighbor(i1);
    iopp = vda->dual_graph_data_structure().mirror_index(f1, i1);
  }
};


//-------------------------------------------------------------------

template<class VDA>
class Find_valid_vertex
{
 public:
  typedef typename VDA::Dual_face_handle  Dual_face_handle;
  typedef std::map<Dual_face_handle,bool> Dual_face_map;

  Dual_face_handle operator()(const VDA* vda,
			      const Dual_face_handle& f) const
  {
    CGAL_precondition( !vda->dual().is_infinite(f) );
    Dual_face_map fmap;
    Dual_face_handle fvalid;
    find_valid_vertex(vda, f, fvalid, fmap);
    CGAL_assertion( fvalid != Dual_face_handle() );
    CGAL_assertion( !vda->dual().is_infinite(fvalid) );
    fmap.clear();
    return fvalid;
  }

 private:
  void find_valid_vertex(const VDA* vda, const Dual_face_handle& cur,
			 Dual_face_handle& fvalid,
			 Dual_face_map& fmap) const
  {
    if ( fmap.find(cur) != fmap.end() ) { return; }
    fmap[cur] = true;


    bool b[3];
    for (int i = 0; i < 3; i++) {
      b[i] = !vda->edge_tester()(cur, i);
    }

    if ( b[0] || b[1] || b[2] ) {
      if ( fvalid == Dual_face_handle() || cur < fvalid ) {
#if 1
	if ( !vda->dual().is_infinite(cur) ) {
	  fvalid = cur;
	}
#else
	fvalid = cur;
#endif
      }
    }

    for (int i = 0; i < 3; i++) {
      if ( !vda->dual().is_infinite(cur->neighbor(i)) && !b[i] ) {
	find_valid_vertex(vda, cur->neighbor(i), fvalid, fmap);
      }
    }
  }

};


//-------------------------------------------------------------------
//-------------------------------------------------------------------

CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_DIAGRAM_2_FINDER_CLASSES_H
