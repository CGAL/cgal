#ifndef CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_GENERIC_H
#define CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_GENERIC_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_triangulation_2.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>

#include <CGAL/draw_triangulation_2.h>
#include <CGAL/internal/Generic_P2T2/Periodic_2_triangulation_vertex_base_2_generic.h>
#include <CGAL/internal/Generic_P2T2/Periodic_2_triangulation_face_base_2_generic.h>
#include <CGAL/utility.h>

#include <utility>
#include <iostream>
#include <fstream>

namespace CGAL {
namespace Periodic_2_triangulations_2 {
namespace internal {

// @todo steal the get_vertex(Face_handle, int, Vertex& /*canonical*/, Offset) from P2T2_GT

// @todo incorporate that in a P2T2GenericTraits class
// @todo number_of_edges/faces + iterators (use an if() to switch between DT2 / P2T2)

// @todo: switch to Bowyer watson?

template < class GT >
class Lattice_2
{
public:
  typedef typename GT::FT                                FT;
  typedef typename GT::Vector_2                          Vector;
  typedef typename GT::Point_2                           Point;

  typedef typename CGAL::Periodic_2_offset_2             Offset;

  typedef std::array<Vector, 2>                        Basis;
  typedef std::array<Vector, 3>                        Voronoi_face_normals;

  // Constructors
  Lattice_2(const Vector& v0, const Vector& v1, const GT& gt = GT())
    : basis_(CGAL::make_array(v0, v1)), gt_(gt)
  {
    initialize();
  }

  Lattice_2(const Basis& basis, const GT& gt = GT())
    : basis_(basis), gt_(gt)
  {
    initialize();
  }

  // Access
  const Basis& basis() const { return basis_; }
  const Voronoi_face_normals& Voronoi_vectors() const { return Vfn_; }
  FT systole_sq_length() const { return systole_sq_length_; }

  // Initialization
  void initialize()
  {
    reduce_basis();
    construct_Voronoi_face_normals();
  }

  void reduce_basis()
  {
    bool reduced = false;
    while (!reduced)
    {
      FT c01 = gt_.compute_scalar_product_2_object()(basis_[0], basis_[1]);
      FT c00 = gt_.compute_scalar_product_2_object()(basis_[0], basis_[0]);
      FT c11 = gt_.compute_scalar_product_2_object()(basis_[1], basis_[1]);
      if (c11 < c00) {
        std::swap(basis_[0], basis_[1]);
        std::swap(c00, c11);
      }
      if (4*c01*c01 <= c00*c00) {
        // Basis is Lagrange-reduced.
        if (c01 > 0) {
          // Negate b1 if necessary to ensure obtuse angle between b0 and b1.
          basis_[1] = gt_.construct_opposite_vector_2_object()(basis_[1]);
        }
        reduced = true;
      } else {
        // Basis is not Lagrange-reduced.
        if (c01 > 0) {
          // b1 -= b0
          basis_[1] = gt_.construct_sum_of_vectors_2_object()(basis_[1],
                        gt_.construct_opposite_vector_2_object()(basis_[0]));
        } else {
          // b1 += b0
          basis_[1] = gt_.construct_sum_of_vectors_2_object()(basis_[1], basis_[0]);
        }
      }
    }
    CGAL_assertion(basis_is_reduced());
  }

  // Only used in assertion check.
  bool basis_is_reduced()
  {
    Vector ext = gt_.construct_opposite_vector_2_object()(
                   gt_.construct_sum_of_vectors_2_object()(basis_[0], basis_[1]));
    return gt_.compute_scalar_product_2_object()(basis_[0], basis_[1]) <= 0 &&
        gt_.compute_scalar_product_2_object()(basis_[0], ext) <= 0 &&
        gt_.compute_scalar_product_2_object()(basis_[1], ext) <= 0;
  }

  void construct_Voronoi_face_normals()
  {
    // @tmp is this really needed or can things be done with predicates?
    Vector third = gt_.construct_opposite_vector_2_object()(
                     gt_.construct_sum_of_vectors_2_object()(basis_[0], basis_[1]));

    Vfn_ = CGAL::make_array(basis_[0], basis_[1], third);

    // @todo check the reduction algorithm, might not be needed to also check the third's length
    systole_sq_length_ = (std::min)((std::min)(basis_[0].squared_length(),
                                               basis_[1].squared_length()),
                                    third.squared_length());
  }

  // Canonicalization
  // @fixme, this shouldn't take the offsetted point, but the canonical point and
  // the offset (to obtain an exact predicate)
  bool is_in_scaled_domain(const Point& p,
                           const FT scaling_factor = 1) const
  {
    for(int i=0; i<3; ++i)
    {
      const Vector& vfn = Vfn_[i];
      const Vector ptv(CGAL::ORIGIN, p);

      const FT sp = gt_.compute_scalar_product_2_object()(ptv, vfn) /
                      gt_.compute_scalar_product_2_object()(vfn, vfn);

      if(!(-0.5 * scaling_factor <= sp && sp < 0.5 * scaling_factor))
        return false;
    }

    return true;
  }

  Point construct_canonical_point(const Point& p) const
  {
    // @check It is fine to do constructions here because an approximation
    // of the exact canonical position of 'p' is fine: we only care about
    // consistency between that approximate position and its offsetted positions

    Point cp = p;

    // @todo factorize it with 'is_in_domain()' (somehow)
    int vfn_pos = 0;
    while(vfn_pos < 3)
    {
      const Vector& vfn = Vfn_[vfn_pos]; // @todo operator(int)
      const Vector ptv(CGAL::ORIGIN, cp);

      const FT sp = gt_.compute_scalar_product_2_object()(ptv, vfn) /
                      gt_.compute_scalar_product_2_object()(vfn, vfn);

      if(-0.5 <= sp && sp < 0.5)
      {
        ++vfn_pos;
      }
      else
      {
        Vector tv = vfn;
        tv = gt_.construct_scaled_vector_2_object()(tv, - std::floor(CGAL::to_double(sp + 0.5) ));
        cp = gt_.construct_translated_point_2_object()(cp, tv);
        vfn_pos = 0;
      }
    }

    return cp;
  }

  Point translate_by_offset(const Point& p, const Offset o) const
  {
    std::cout << "Reduced basis b[0] = " << basis_[0] << std::endl;
    std::cout << "Reduced basis b[1] = " << basis_[1] << std::endl;

    std::cout << "translate_by_offset(" << p << " Off: " << o << ") = ";

    Vector translation = gt_.construct_sum_of_vectors_2_object()(
        gt_.construct_scaled_vector_2_object()(basis_[0], o.x()),
        gt_.construct_scaled_vector_2_object()(basis_[1], o.y()));

    std::cout << gt_.construct_translated_point_2_object()(p, translation) << std::endl;

    return gt_.construct_translated_point_2_object()(p, translation);
  }

private:
  FT systole_sq_length_;
  std::array<Vector, 2> basis_;
  std::array<Vector, 3> Vfn_;
  const GT gt_; // @todo pointer
};

template < typename K_, typename Construct_point_2_base_>
class Lattice_construct_point_2
  : public Construct_point_2_base_
{
  typedef Construct_point_2_base_            Base;
  typedef K_                                 Kernel;

  typedef typename Kernel::Point_2           Point;
  typedef CGAL::Periodic_2_offset_2          Offset;

  typedef internal::Lattice_2<K_>            Lattice;

public:
  Lattice_construct_point_2(const Lattice* lattice, const Base& cp)
    : Base(cp), lattice_(lattice)
  { }

  using Base::operator();

  Point operator()(const Point& p, const Offset& o) const
  {
    CGAL_assertion(lattice_ != nullptr);
    return lattice_->translate_by_offset(p, o);
  }

private:
  const Lattice* lattice_;
};

} // end namespace internal
} // end namespace Periodic_2_triangulations_2

template < class GT,
           class TDS = Triangulation_data_structure_2 <
                         Periodic_2_triangulation_vertex_base_2_generic<GT>,
                         Periodic_2_triangulation_face_base_2_generic<GT> > >
class Periodic_2_Delaunay_triangulation_2_generic
{
  typedef Periodic_2_Delaunay_triangulation_2_generic<GT, TDS>          Self;

public:
  typedef TDS                                  Triangulation_data_structure;
  typedef GT                                   Geom_traits;

  typedef typename CGAL::Periodic_2_offset_2   Offset; // @fixme should be defined by the traits

  typedef typename GT::FT                      FT;
  typedef typename GT::Point_2                 Point;
  typedef typename GT::Segment_2               Segment;
  typedef typename GT::Vector_2                Vector;
  typedef typename GT::Triangle_2              Triangle;

  typedef std::pair<Point, Offset>              Periodic_point;
  typedef array< std::pair<Point, Offset>, 2>   Periodic_segment;
  typedef array< std::pair<Point, Offset>, 3>   Periodic_triangle;
  typedef array< std::pair<Point, Offset>, 4>   Periodic_tetrahedron;

  typedef typename TDS::size_type               size_type;

  typedef typename TDS::Face_handle             Face_handle;
  typedef typename TDS::Vertex_handle           Vertex_handle;
  typedef typename TDS::Edge                    Edge;

  typedef typename TDS::Vertex_circulator       Vertex_circulator;
  typedef typename TDS::Face_circulator         Face_circulator;

  //Tag to distinguish Delaunay from regular triangulations
  typedef Tag_false                             Weighted_tag;

  // Tag to distinguish periodic triangulations from others
  typedef Tag_true                              Periodic_tag;

  typedef Periodic_2_triangulations_2::internal::Lattice_2<GT> Lattice;
  typedef typename Lattice::Basis                              Basis;
  typedef typename Lattice::Voronoi_face_normals               Voronoi_face_normals;

  typedef CGAL::Delaunay_triangulation_2<GT, TDS>       DT2;

  typedef typename GT::Construct_point_2                         Base_CP2;
  typedef Periodic_2_triangulations_2::internal::Lattice_construct_point_2<GT, Base_CP2> CP2;

  typedef CGAL::Periodic_2_triangulation_traits_base_2<GT, Offset, Lattice, CP2> P2T2_GT;
  typedef CGAL::Periodic_2_triangulation_2<P2T2_GT, TDS>                    P2T2;

  enum Locate_type
  {
    VERTEX = 0,
    EDGE, //1
    FACE, //2
    EMPTY , //4
    OUTSIDE_CONVEX_HULL, // unused, for compatibility with Alpha shapes
    OUTSIDE_AFFINE_HULL // unused, for compatibility with Alpha shapes
  };

  /// Constructors
  template <class InputIterator>
  Periodic_2_Delaunay_triangulation_2_generic(InputIterator first, InputIterator beyond,
                                              const Basis& basis,
                                              const GT& gt = GT())
    : lattice_(basis), is_1_cover(false), gt_(gt), p2t2_gt_(lattice_, gt_), p2t2(lattice_, p2t2_gt_)
  {
    sq_circumradius_threshold = 0.25 * lattice_.systole_sq_length();
    insert(first, beyond);
  }

#ifndef CGAL_CFG_NO_CPP0X_DELETED_AND_DEFAULT_FUNCTIONS
  Periodic_2_Delaunay_triangulation_2_generic& operator=(const Periodic_2_Delaunay_triangulation_2_generic&)=default;
#endif

  const Geom_traits& geom_traits() const { return gt_; }
  const Triangulation_data_structure & tds() const { return tds_; }

// @todo
//  void copy_triangulation(const Periodic_2_Delaunay_triangulation_2_generic &tr) { }
//  void swap(Periodic_2_Delaunay_triangulation_2_generic &tr) { }
//  void clear() { }

  // number of canonical simplicies of each dimension

  /// Returns the number of vertices. Counts all vertices that are
  /// representatives of the same point in the 1-cover as one vertex.
  size_type number_of_vertices() const
  {
    if(is_1_cover)
      return p2t2.number_of_vertices();
    else
      return dt2.number_of_vertices() / 9;
  }

  size_type number_of_edges() const
  {
    if(is_1_cover)
      return p2t2.number_of_edges();
    else
    {
      // Exploiting Euler's formula that #f - #e + #v = 0 on the torus
      return number_of_faces() + number_of_vertices();

    //   // Alternative naive implementation
    //   size_type ne = 0;
    //   typename DT2::Finite_edges_iterator eit = dt2.finite_edges_begin(),
    //                                       eend = dt2.finite_edges_end();
    //   for(; eit!=eend; ++eit)
    //     if(is_canonical(*eit))
    //       ++ne;

    //   return ne;
    }
  }

  size_type number_of_faces() const
  {
    if(is_1_cover)
      return p2t2.number_of_faces();
    else
    {
      // @todo something better than this naive way (flag edge/face classes)
      size_type nf = 0;
      typename DT2::Finite_faces_iterator fit = dt2.finite_faces_begin(),
                                          fend = dt2.finite_faces_end();
      for(; fit!=fend; ++fit)
        if(is_canonical(fit))
          ++nf;

      return nf;
    }
  }

//   Vertices_iterator vertices_begin() const { }
//   Vertices_iterator vertices_end() const { }
//   Edges_iterator edges_begin() const { }
//   Edges_iterator edges_end() const { }
//   Faces_iterator faces_begin() const { }
//   Faces_iterator faces_end() const { }

  // (and finite versions, e.g. Finite_vertices_iterator finite_vertices_iterator, for other packages)

//   Vertex_circulator adjacent_vertices(Vertex_handle vh) const { }
//   Edge_circulator incident_edges(Vertex_handle vh) const { }
//   Face_circulator incident_faces(Vertex_handle vh) const { }

  int dimension() const { return (number_of_vertices() == 0) ? -2 : 2; }

  // @todo two versions, one using the sufficient condition, one checking for real.
  // "For real" --> P3T3 is_embeddable_in_...
  // "For real" is implemented below.
  bool is_simplicial_complex() const
  {
    // Ensure there is no edge between a vertex and itself.
    typename DT2::Finite_faces_iterator fit = dt2.faces_begin(),
                                 fend = dt2.faces_end();
    for(; fit!=fend; ++fit) {
      if (!is_canonical(fit))
        continue;
      Vertex_handle vh0 = canonical_vertex(fit->vertex(0));
      Vertex_handle vh1 = canonical_vertex(fit->vertex(1));
      Vertex_handle vh2 = canonical_vertex(fit->vertex(2));

      if (vh0 == vh1 || vh0 == vh2 || vh1 == vh2) {
        return false;
      }
    }

    // Ensure there are no two edges between the same pair of vertices.
    typename DT2::Finite_vertices_iterator vit = dt2.vertices_begin(),
                                 vend = dt2.vertices_end();
    for(; vit!=vend; ++vit) {
      if (!is_canonical(vit))
        continue;

      typename DT2::Vertex_circulator vc = dt2.incident_vertices(vit), done = vc;
      std::unordered_set<Vertex_handle> neighbours;
      do
      {
        Vertex_handle cv = canonical_vertex(vc);
        if (neighbours.find(cv) != neighbours.end())
          return false; // Some neighbouring vertex appeared multiple times.
        neighbours.insert(cv);
      }
      while(++vc != done);
    }

    return true;
  }

  /// Constructions
  Point construct_point(const Point& /*p*/, const Offset& /*off*/) const { }

  // @todo is this useful?
  std::pair<Point /*canonical point*/, Offset> periodic_point(const Point& /*p*/) const { }

  /// Canonicalization
  Point construct_canonical_point(const Point& p) const
  {
    return lattice_.construct_canonical_point(p);
  }

  template <class InputIterator>
  std::vector<Point> construct_canonical_points(InputIterator first, InputIterator beyond) const
  {
    std::vector<Point> canonical_points;

    while(first != beyond)
    {
      const Point& p = *first++;
      canonical_points.push_back(construct_canonical_point(p));
    }

    return canonical_points;
  }

  FT compute_squared_circumradius(const Face_handle fh) const
  {
    return geom_traits().compute_squared_radius_2_object()(dt2.point(fh, 0),
                                                           dt2.point(fh, 1),
                                                           dt2.point(fh, 2));
  }

  /// Canonicity functions
  bool is_canonical(const Point& p) const
  {
    return lattice_.is_in_scaled_domain(p, 1);
  }

  Vertex_handle canonical_vertex(const Vertex_handle vh) const
  {
    return (is_canonical(vh) ? vh : canonical_vertices.at(vh));
  }

  bool is_canonical(const Vertex_handle vh) const
  {
    return (vh->offset() == Offset(0, 0));
  }

  Offset compute_offset(const Edge e) const
  {
    Face_handle fh = e.first;
    int i = e.second;
    return compute_offset(fh->vertex(dt2.cw(i)), fh->vertex(dt2.ccw(i)));
  }

  bool is_canonical(const Edge e) const
  {
    if(dt2.is_infinite(e.first))
      return false;

    return compute_offset(e) == Offset(0, 0);
  }

  // Offset of an edge represented by its two vertices
  Offset compute_offset(const Vertex_handle vh1, const Vertex_handle vh2) const
  {
    return min(vh1->offset(), vh2->offset());
  }

  // Canonicity of an edge represented by its two vertices
  bool is_canonical(const Vertex_handle vh1, const Vertex_handle vh2) const
  {
    if(dt2.is_infinite(vh1) || dt2.is_infinite(vh2))
      return false;

    return compute_offset(vh1, vh2) == Offset(0, 0);
  }

  Offset compute_offset(const Face_handle fh) const
  {
    return min(fh->vertex(0)->offset(), fh->vertex(1)->offset(), fh->vertex(2)->offset());
  }

  bool is_canonical(const Face_handle fh) const
  {
    if(dt2.is_infinite(fh))
      return false;

    return compute_offset(fh) == Offset(0, 0);
  }

  template <typename ForwardFaceIterator>
  void mark_canonical_faces(ForwardFaceIterator fit, ForwardFaceIterator beyond)
  {
    while(fit != beyond)
    {
      fit->set_canonical_flag(is_canonical(fit++));
    }
  }

  void reset_all_canonicity()
  {
    typename DT2::Finite_faces_iterator fit = dt2.faces_begin(),
                                 fend = dt2.faces_end();
    for(; fit!=fend; ++fit)
      fit->set_canonical_flag(false);
  }

  void mark_canonical_faces()
  {
    return mark_canonical_faces(dt2.finite_faces_begin(), dt2.finite_faces_end());
  }

  void mark_canonical_faces(Vertex_handle vh)
  {
    CGAL_precondition(vh != Vertex_handle());

    if(dt2.dimension() != 2)
      return;

    typename DT2::Face_circulator fc = dt2.incident_faces(vh), done = fc;
    do
    {
      fc->set_canonical_flag(is_canonical(fc));
    }
    while(++fc != done);
  }

  Face_handle get_canonical_face(const Face_handle& fh) const
  {
    std::cout << "Getting canonical face of: " << dt2.point(fh, 0) << " "
                                               << dt2.point(fh, 1) << " "
                                               << dt2.point(fh, 2) << std::endl;

    return find_translated_face(fh, -compute_offset(fh));
  }

  /// Low level functions to mascarade the DT2 as a periodic triangulation
  Point construct_barycenter(const Face_handle& fh) const
  {
    Vector v0 = Vector(CGAL::ORIGIN, fh->vertex(0)->point());
    Vector v1 = Vector(CGAL::ORIGIN, fh->vertex(1)->point());
    Vector v2 = Vector(CGAL::ORIGIN, fh->vertex(2)->point());
    Vector bcv = gt_.construct_scaled_vector_2_object()(
              gt_.construct_sum_of_vectors_2_object()(v0,
                gt_.construct_sum_of_vectors_2_object()(v1,v2)),
            FT(1)/3);
    return Point(bcv.x(), bcv.y());
  }

  // Given a face having a vertex in the domain, and an offset such that
  // at least one of the translated vertices has a vertex in the domain,
  // return the handle of the translated face.
  Face_handle find_translated_face(const Face_handle& fh, const Offset& o) const
  {
    // The code commented out below does the same and is simpler,
    // but uses a construction and point location.
//    Point translated_barycenter = lattice_.translate_by_offset(
//        construct_barycenter(fh), o);
//    return dt2.locate(translated_barycenter);

    // Find a vertex whose translate is in the domain.
    bool vertex_found = false;
    int j=0;
    for (; j<3; ++j) {
      if (fh->vertex(j)->offset() == -o) {
        vertex_found = true;
        break;
      }
    }
    CGAL_assertion(vertex_found);
    Vertex_handle cv = canonical_vertices.at(fh->vertex(j));
    Vertex_handle v_ccw = fh->vertex(dt2.ccw(j));
    Vertex_handle v_cw = fh->vertex(dt2.cw(j));

    // Scan through the incident faces and find the one that is
    // equivalent to fh.
    Face_circulator fc = dt2.incident_faces(cv), done = fc;
    do
    {
      if(dt2.is_infinite(fc)) // shouldn't ever happen
        continue;

      int cj = fc->index(cv);
      Vertex_handle cv_ccw = fc->vertex(dt2.ccw(cj));
      Vertex_handle cv_cw = fc->vertex(dt2.cw(cj));

      if (canonical_vertices.at(cv_ccw) == canonical_vertices.at(v_ccw)
          && cv_ccw->offset() == v_ccw->offset() + o
          && canonical_vertices.at(cv_cw) == canonical_vertices.at(v_cw)
          && cv_cw->offset() == v_cw->offset() + o) {
        return Face_handle(fc);
      }
    }
    while(++fc != done);
    CGAL_assertion_msg(false, "couldn't find face");
    return Face_handle();
  }

  // @todo something smarter, this function is core to everything else
  // How could we keep "canonical neighbors" in memory instead of having
  // to find them later...?
  Face_handle neighbor(Face_handle fh, int i) const
  {
    if(is_canonical(fh->neighbor(i)))
      return fh->neighbor(i);

    // Translate the face so that the corresponding edge is canonical.
    Offset edge_off = compute_offset(std::make_pair(fh, i));
    Face_handle fh_trans = find_translated_face(fh, -edge_off);

    // Find the vertex corresponding to vertex(i) in the original.
    bool vertex_found = false;
    int j=0;
    for (; j<3; ++j) {
      if (canonical_vertices.at(fh_trans->vertex(j)) == canonical_vertices.at(fh->vertex(i))
          && fh_trans->vertex(j)->offset() == fh->vertex(i)->offset() - edge_off) {
        vertex_found = true;
        break;
      }
    }
    CGAL_assertion(vertex_found);

    // Get the neighbour in DT2 and check if it's canonical.
    Face_handle neighbor = fh_trans->neighbor(j);
    if(is_canonical(neighbor))
      return neighbor;
    else // if not, find the canonical translate.
      return get_canonical_face(neighbor);
  }

  /// Iterators and Circulators
  // @todo this should return a (custom) circulator
  std::set<Face_handle> incident_faces(Vertex_handle vh) const
  {
    std::set<Face_handle> ifhs;

    // Need to start from the canonical vertex to make sure there is no "weird" (to be defined properly) faces
    if(!is_canonical(vh))
      vh = canonical_vertices.at(vh);

    Face_circulator tds_fc = dt2.incident_faces(vh), done = tds_fc;

    do
    {
      ifhs.insert(get_canonical_face(tds_fc));

      // @todo when proper periodic traits are used to construct dt2 (that is,
      // insert (p, off) instead of constructing the offsetted point,
      // then the vertices will store the _canonical_ (geometric) point,
      // and we don't have to use the canonical_vertices[] but rather we can
      // just compare fh->vertex(0)->point() and vh->point()
      //
      // "Real" offsetted point is then done via a 'construct_point(cp, off)' function
      // in the triangulation class

      // The alternative is to use the neighbourhood relation of the faces, and
      // traverse along edges around the central vertex, in the process possibly
      // shifting the center vertex around to stay canonical.
      // Incomplete code of this in a previous commit.
    }
    while(++tds_fc != done);

    return ifhs;
  }

  /// Locate functions
//  Face_handle locate(const Point& p, Offset& lo,
//                     Locate_type& lt, int& li,
//                     Face_handle start = Face_handle()) const
//  {

//  }

//  Face_handle locate(const Point& p,
//                     Locate_type& lt, int& li,
//                     Face_handle start = Face_handle()) const
//  {
//    Offset lo;
//    return locate(p, lo, lt, li, start);
//  }

//  Face_handle locate(const Point& p,
//                     Face_handle start = Face_handle()) const
//  {
//    Locate_type lt;
//    int li;
//    return locate(p, lt, li, start);
//  }

  /// Insertion and removal
  Vertex_handle insert(const Point& p)
  {
    if(is_1_cover)
      return insert_in_p2t2(p);
    else
      return insert_in_dt2(p);
  }

  Vertex_handle insert_in_dt2(const Point& p)
  {
    const Point cp = lattice_.construct_canonical_point(p);
    std::cout << "Insert (DT2): " << p << " canonical: " << cp << std::endl;

    std::cout << dt2.number_of_vertices() << " vertices" << std::endl;
    if(dt2.dimension() >= 2) // equivalent to !dt2.empty() since we insert duplicate vertices
    {
      // @todo avoid recomputing the conflict zone if possible (done also 'insert', sort of)
      std::vector<Face_handle> faces_in_conflict;
      dt2.get_conflicts(cp, std::back_inserter(faces_in_conflict));
      std::cout << faces_in_conflict.size() << " faces in conflict" << std::endl;

      for(Face_handle fh : faces_in_conflict)
      {
        // @fixme? is this safe? Are all faces returned by the get_conflicts
        // faces of the periodic triangulation and not "boundary" faces
        faces_with_too_big_circumradius.erase(get_canonical_face(fh));
      }
    }

    Vertex_handle vh = dt2.insert(cp);

    CGAL_assertion(vh != Vertex_handle());
    vh->set_offset(Offset(0, 0));

    mark_canonical_faces(vh);

    canonical_vertices[vh] = vh;
    for (const std::vector<int> off : overlapping_offsets)
    {
      // @fixme Insert without constructions (cf P3T3)
      const Vector off_v = gt_.construct_sum_of_vectors_2_object()(
                             gt_.construct_scaled_vector_2_object()(
                               lattice_.basis()[0], off[0]),
                             gt_.construct_scaled_vector_2_object()(
                               lattice_.basis()[1], off[1]));

      const Point off_p = cp + off_v;

      if(lattice_.is_in_scaled_domain(off_p, 3))
      {
        Vertex_handle vh_copy = dt2.insert(off_p);
        CGAL_assertion(vh_copy != Vertex_handle());
        vh_copy->set_offset(Offset(off[0], off[1]));

        canonical_vertices[vh_copy] = vh;

        mark_canonical_faces(vh_copy);
      }
    }

    // Update the current maximum circumradius value
    std::cout << "Gather faces with too big circumradius" << std::endl;
    Face_circulator fc = dt2.incident_faces(vh), done(fc);
    do
    {
      CGAL_assertion(!dt2.is_infinite(fc));

      Face_handle cfh = get_canonical_face(fc);
      const FT sq_cr = compute_squared_circumradius(cfh);
      std::cout << " sq_cr: " << sq_cr << " sys:" << sq_circumradius_threshold << std::endl;
      if(sq_cr > sq_circumradius_threshold)
        faces_with_too_big_circumradius.insert(cfh);
    }
    while(++fc != done);

    std::cout << faces_with_too_big_circumradius.size() << " faces with too big sq_cr" << std::endl;
    // if(faces_with_too_big_circumradius.empty())
    //   convert_to_1_cover();

    return vh;
  }

  Vertex_handle insert_in_p2t2(const Point& p)
  {
    const Point cp = lattice_.construct_canonical_point(p);
    std::cout << "Insert (P2T2): " << p << " canonical: " << cp << std::endl;

    return p2t2.insert(p);
  }

  template <class InputIterator>
  void insert(InputIterator first, InputIterator beyond)
  {
    // @todo mirror it with other insert functions (e.g. spatial sorting)

    std::size_t np = std::distance(first, beyond);

    while(first != beyond)
      insert(*first++);

    CGAL_postcondition(dt2.dimension() == 2);
    CGAL_postcondition(dt2.number_of_vertices() == 9 * np);
    CGAL_postcondition(dt2.is_valid());
  }

  void add_edge_to_incident_faces_map(const Face_handle fh, int i,
                                      std::map<std::set<Vertex_handle>,
                                               std::vector<std::pair<Face_handle, int> > >& incident_faces_map)
  {
    typedef std::set<Vertex_handle>                                            Edge_vertices;
    typedef std::pair<Face_handle, int>                                        Incident_face;
    typedef std::map<Edge_vertices, std::vector<Incident_face> >               Incident_faces_map;

    // the opposite vertex of f in c is i
    Edge_vertices e;
    e.insert(fh->vertex((i + 1) % 3));
    e.insert(fh->vertex((i + 2) % 3));
    CGAL_precondition(e.size() == 2);

    Incident_face icf = std::make_pair(fh, i);
    std::vector<Incident_face> vec;
    vec.push_back(icf);

    std::cout << "set incidences for edge: " << std::endl;
    for(const auto v : e)
      std::cout << v->point() << " ";
    std::cout << std::endl;

    std::pair<typename Incident_faces_map::iterator, bool> is_insert_successful =
        incident_faces_map.insert(std::make_pair(e, vec));
    if(!is_insert_successful.second) // the entry already exists in the map
    {
      // a facet must have exactly two incident faces
      CGAL_assertion(is_insert_successful.first->second.size() == 1);
      is_insert_successful.first->second.push_back(icf);
    }
  }

  void convert_to_1_cover()
  {
    is_1_cover = true;

    p2t2.clear();
    p2t2.tds().set_dimension(2);

    std::unordered_map<Vertex_handle /*dt2*/, Vertex_handle /*p2t2*/> vertex_correspondence_map;

    typename DT2::Finite_vertices_iterator vit = dt2.finite_vertices_begin(),
                                           vend = dt2.finite_vertices_end();
    for(; vit!=vend; ++vit)
    {
      if(!is_canonical(vit))
        continue;

      Vertex_handle vh = p2t2.tds().create_vertex();
      vh->set_point(vit->point());
      vertex_correspondence_map[vit] = vh;
    }

    // @todo array instead of vector
    typedef std::map<std::set<Vertex_handle>, // two vertices of an edge
                     std::vector<std::pair<Face_handle, int> > > Incident_faces_map;
    Incident_faces_map incident_faces_map;

    size_type cfc = 0;
    typename DT2::Finite_faces_iterator fit = dt2.finite_faces_begin(),
                                        fend = dt2.finite_faces_end();
    for(; fit!=fend; ++fit)
    {
      if(!is_canonical(fit))
        continue;

      ++cfc;

      Vertex_handle p2t2_vh0 = vertex_correspondence_map.at(canonical_vertex(fit->vertex(0)));
      Vertex_handle p2t2_vh1 = vertex_correspondence_map.at(canonical_vertex(fit->vertex(1)));
      Vertex_handle p2t2_vh2 = vertex_correspondence_map.at(canonical_vertex(fit->vertex(2)));

      Face_handle fh = p2t2.tds().create_face(p2t2_vh0, p2t2_vh1, p2t2_vh2);

      fh->set_offsets(fit->vertex(0)->offset()-fit->vertex(0)->offset(),
                      fit->vertex(1)->offset()-fit->vertex(0)->offset(),
                      fit->vertex(2)->offset()-fit->vertex(0)->offset());

      add_edge_to_incident_faces_map(fh, 0, incident_faces_map);
      add_edge_to_incident_faces_map(fh, 1, incident_faces_map);
      add_edge_to_incident_faces_map(fh, 2, incident_faces_map);

      // Set up incident face information
      for(int i=0; i<3; ++i)
      {
        if(fh->vertex(i)->face() == Face_handle())
          fh->vertex(i)->set_face(fh);
      }
    }

    std::cout << dt2.number_of_vertices() / 9 << " canonical vertices in dt2" << std::endl;
    std::cout << cfc << " canonical faces in dt2" << std::endl;

    // Set up adjacencies
    typename Incident_faces_map::const_iterator ifit = incident_faces_map.begin();
    for(; ifit!=incident_faces_map.end(); ++ifit)
    {
      const std::vector<std::pair<Face_handle, int> >& adjacent_faces = ifit->second;
      CGAL_assertion(adjacent_faces.size() == 2);

      Face_handle f0 = adjacent_faces[0].first;
      int i0 = adjacent_faces[0].second;
      Face_handle f1 = adjacent_faces[1].first;
      int i1 = adjacent_faces[1].second;

      p2t2.tds().set_adjacency(f0, i0, f1, i1);
    }

    std::cout << p2t2.number_of_vertices() << " vertices in p2t2" << std::endl;
    std::cout << p2t2.number_of_faces() << " canonical faces in p2t2" << std::endl;

    CGAL_postcondition(p2t2.tds().is_valid(true));
    CGAL_postcondition(p2t2.is_valid(true));
  }

  void draw_p2t2() const
  {
    std::ofstream out("p2t2.off");
    out << "OFF\n";

    std::map<Point, int> ids;
    std::vector<std::array<int, 3> > faces;

    int pid = 0;
    for(auto fit=p2t2.faces_begin(); fit!=p2t2.faces_end(); ++fit)
    {
      std::array<int, 3> face;
      for(int i=0; i<3; ++i)
      {
        auto itb = ids.insert(std::make_pair(p2t2.point(fit, i), pid));
        if(itb.second)
          ++pid;
        face[i] = itb.first->second;
      }
      faces.push_back(face);
    }

    CGAL_assertion(faces.size() == p2t2.number_of_faces());

    out << ids.size() << " " << p2t2.number_of_faces() << " 0\n";

    for(const auto& e : ids)
      out << e.first << " 0" << std::endl;

    for(const auto& face : faces)
      out << "3 " << face[0] << " " << face[1] << " " << face[2] << std::endl;

    out << std::endl;
  }

public:
  Geom_traits gt_;
  DT2 dt2; // @tmp, shouldn't be exposed
  P2T2_GT p2t2_gt_;
  P2T2 p2t2;

private:
  Lattice lattice_; // @tmp to be incorporated in the geometric traits

  bool is_1_cover;

  std::unordered_map<Vertex_handle /*periodic copy*/,
                     Vertex_handle /*canonical*/> canonical_vertices;

  std::set<Face_handle> faces_with_too_big_circumradius;
  FT sq_circumradius_threshold;

  Triangulation_data_structure tds_;

  // A list of those offsets such that the domain translated along the offset
  // overlaps the scaled domain.
  // Could be static?
  const std::vector<std::vector<int> > overlapping_offsets = {
      // entirely contained in scaled domains
      {-1, -1}, {0, 1}, {1, 0}, {-1, 0}, {0, -1}, {1, 1},
      // intersecting the scaled domain
      {-1, -2}, {1, 2}, {-2, -1}, {2, 1}, {-1, 1}, {1, -1}
    };
};

} //namespace CGAL

#endif // CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_GENERIC_H
