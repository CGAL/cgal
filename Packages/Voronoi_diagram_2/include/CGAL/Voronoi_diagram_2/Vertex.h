#ifndef CGAL_VORONOI_DIAGRAM_2_VERTEX_H
#define CGAL_VORONOI_DIAGRAM_2_VERTEX_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Finder_classes.h>


CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

template<class VDA>
class Vertex
{
 private:
  typedef Vertex<VDA>                     Self;
  typedef Triangulation_cw_ccw_2          CW_CCW_2;

  typedef typename VDA::Dual_face_handle  Dual_face_handle;

 private:
  Dual_face_handle find_valid_vertex(const Dual_face_handle& f) const
  {
    return Find_valid_vertex<VDA>()(vda_, f);
  }

 public:
  typedef typename VDA::Halfedge          Halfedge;
  typedef typename VDA::Halfedge_handle   Halfedge_handle;
  typedef typename VDA::Face              Face;
  typedef typename VDA::Face_handle       Face_handle;
  typedef typename VDA::Vertex_handle     Vertex_handle;
  typedef typename VDA::size_type         size_type;

  typedef typename VDA::Halfedge_around_vertex_circulator
  Halfedge_around_vertex_circulator;

 public:
  Vertex(const VDA* vda = NULL) : vda_(vda) {}
  Vertex(const VDA* vda, Dual_face_handle f) : vda_(vda), f_(f) {
    CGAL_precondition( !vda_->dual().is_infinite(f_) );
  }

  Halfedge_handle halfedge() const {
    Dual_face_handle fvalid = find_valid_vertex(f_);
    for (int i = 0; i < 3; i++) {
      int ccw_i = CW_CCW_2::ccw(i);
      
      // if I want to return also infinite edges replace the test in
      // the if statement by the following test (i.e., should omit the
      // testing for infinity):
      //           !vda_->edge_tester()(vda_->dual(), fvalid, i)
      if ( !vda_->edge_tester()(vda_->dual(), fvalid, i) &&
	   !vda_->dual().is_infinite(fvalid, i) ) {
	if ( vda_->face_tester()(vda_->dual(), vda_->edge_tester(),
				 fvalid->vertex(ccw_i)) ) {
	  Dual_face_handle fopp;
	  int iopp, i_mirror = vda_->dual().tds().mirror_index(fvalid, i);

	  Find_opposite_halfedge<VDA>()(vda_,
					fvalid->neighbor(i),
					i_mirror,
					fopp, iopp);
#ifndef CGAL_NO_ASSERTIONS
	  Halfedge h(vda_, fopp, iopp);
	  Vertex_handle v_this(*this);
	  CGAL_assertion( h.target() == v_this );
#endif
	  return Halfedge_handle( Halfedge(vda_, fopp, iopp) );
	} else {
#ifndef CGAL_NO_ASSERTIONS
	  Halfedge h(vda_, fvalid, i);
	  Vertex_handle v_this(*this);
	  CGAL_assertion( h.target() == v_this );
#endif
	  return Halfedge_handle( Halfedge(vda_, fvalid, i) );
	}
      }
    }

    bool this_line_should_never_have_been_reached = false;
    CGAL_assertion(this_line_should_never_have_been_reached);
    return Halfedge_handle();
  }

  Halfedge_around_vertex_circulator incident_halfedges() const {
    CGAL_assertion( halfedge()->target() == Vertex_handle(*this) );
    return Halfedge_around_vertex_circulator( halfedge() );
  }

  bool is_incident_edge(const Halfedge_handle& he) const {
    return he->target() == Vertex_handle(*this) ||
      he->source() == Vertex_handle(*this);
  }

  bool is_incident_face(const Face_handle& f) const {
    Halfedge_around_vertex_circulator hc = incident_halfedges();
    Halfedge_around_vertex_circulator hc_start = hc;
    do {
      if ( (*hc)->face() == f ) { return true; }
      ++hc;
    } while ( hc != hc_start );
    return false;
  }

  size_type degree() const {
    Halfedge_around_vertex_circulator hc = incident_halfedges();
    Halfedge_around_vertex_circulator hc_start = hc;
    size_type deg = 0;
    do {
      hc++;
      deg++;
    } while ( hc != hc_start );
    return deg;
  }

  bool operator==(const Self& other) const {
    if ( vda_ == NULL ) { return other.vda_ == NULL; }
    if ( other.vda_ == NULL ) { return vda_ == NULL; }
    return ( vda_ == other.vda_ && f_ == other.f_ );
  }

  bool operator!=(const Self& other) const {
    return !((*this) == other);
  }

  // temporary
  const Dual_face_handle& dual_face() const { return f_; }

  bool is_valid() const {
    if ( vda_ == NULL ) { return true; }

    bool valid = !vda_->dual().is_infinite(f_);
    // THE FOLLOWING LINE CREATES A PROBLEM FOR SOME REASON...
    //    valid = valid && is_incident_edge( halfedge() );

    Vertex_handle v_this(*this);

    valid = valid && halfedge()->target() == v_this;
    valid = valid && halfedge()->opposite()->source() == v_this;

    Halfedge_around_vertex_circulator hc = incident_halfedges();
    Halfedge_around_vertex_circulator hc_start = hc;
    do {
      valid = valid && (*hc)->target() == v_this;
      ++hc;
    } while ( hc != hc_start );

    Halfedge_handle hhc = *incident_halfedges();
    Halfedge_handle hhc_start = hhc;
    do {
      valid = valid && hhc->target() == v_this;
      hhc = hhc->next()->twin();
    } while ( hhc != hhc_start );

    return valid;
  }

 private:
  const VDA* vda_;
  Dual_face_handle f_;
};


CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE


#endif // CGAL_VORONOI_DIAGRAM_2_VERTEX_H
