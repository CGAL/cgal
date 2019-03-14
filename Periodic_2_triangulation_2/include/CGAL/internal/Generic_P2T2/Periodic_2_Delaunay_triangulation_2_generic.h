#ifndef CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_GENERIC_H
#define CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_GENERIC_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>

#include <CGAL/internal/Generic_P2T2/Periodic_2_triangulation_vertex_base_2_generic.h>
#include <CGAL/internal/Generic_P2T2/Periodic_2_triangulation_face_base_2_generic.h>
#include <CGAL/utility.h>

#include <CGAL/draw_triangulation_2.h>
#include <CGAL/unordered.h>

#include <utility>
#include <iostream>

namespace CGAL {

namespace Periodic_2_triangulations_2 {

namespace internal {

// @todo incorporate that in a P2T2GenericTraits class
template < class GT >
class Lattice_2
{
public:
  typedef typename GT::FT                                FT;
  typedef typename GT::Vector_2                          Vector;
  typedef typename GT::Point_2                           Point;

  typedef cpp11::array<Vector, 2>                        Basis;
  typedef cpp11::array<Vector, 3>                        Voronoi_face_normals;

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

  // Initialization
  void initialize()
  {
    reduce_basis();
    construct_Voronoi_face_normals();
  }

  // @tmp
  void reduce_basis()
  {
  }

  void construct_Voronoi_face_normals()
  {
    // @tmp
    Vector third = gt_.construct_opposite_vector_2_object()(
                     gt_.construct_sum_of_vectors_2_object()(basis_[0], basis_[1]));

    Vfn_ = CGAL::make_array(basis_[0], basis_[1], third);
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
        tv = gt_.construct_scaled_vector_2_object()(tv, - std::floor(sp + 0.5) );
        cp = gt_.construct_translated_point_2_object()(cp, tv);
        vfn_pos = 0;
      }
    }

    return cp;
  }

private:
  CGAL::cpp11::array<Vector, 2> basis_;
  CGAL::cpp11::array<Vector, 3> Vfn_;
  const GT& gt_;
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

  typedef typename CGAL::Periodic_2_offset_2   Offset;

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
  typedef typename Lattice::Basis                             Basis;
  typedef typename Lattice::Voronoi_face_normals              Voronoi_face_normals;

  typedef CGAL::Delaunay_triangulation_2<GT, TDS>       DT2;

  typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<GT, Offset> P2T2_GT;
  typedef CGAL::Periodic_2_triangulation_2<P2T2_GT, TDS>               P2T2;

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
    : lattice_(basis), is_1_cover(false), gt_(gt)
  {
    insert(first, beyond);
  }

#ifndef CGAL_CFG_NO_CPP0X_DELETED_AND_DEFAULT_FUNCTIONS
  Periodic_2_Delaunay_triangulation_2_generic& operator=(const Periodic_2_Delaunay_triangulation_2_generic&)=default;
#endif

  const Geom_traits& geom_traits() const { return gt_; }
  const Triangulation_data_structure & tds() const { return tds_; }

//  void copy_triangulation(const Periodic_2_Delaunay_triangulation_2_generic &tr) { }
//  void swap(Periodic_2_Delaunay_triangulation_2_generic &tr) { }
//  void clear() { }

  void add_edge_to_incident_faces_map(const Face_handle fh, int i,
                                      boost::unordered_map<
                                        std::set<Vertex_handle>,
                                        std::vector< // @todo array
                                          std::pair<Face_handle, int> > >& incident_faces_map)
  {
    typedef std::set<Vertex_handle>                                            Edge_vertices;
    typedef std::pair<Face_handle, int>                                        Incident_face;
    typedef boost::unordered_map<Edge_vertices, std::vector<Incident_face> >   Incident_faces_map;

    // the opposite vertex of f in c is i
    Edge_vertices e;
    e.insert(fh->vertex((i + 1) % 3));
    e.insert(fh->vertex((i + 2) % 3));
    CGAL_precondition(e.size() == 2);

    Incident_face icf = std::make_pair(fh, i);
    std::vector<Incident_face> vec;
    vec.push_back(icf);

    std::pair<typename Incident_faces_map::iterator, bool> is_insert_successful =
        incident_faces_map.insert(std::make_pair(e, vec));
    if(!is_insert_successful.second) // the entry already exists in the map
    {
      // a facet must have exactly two incident faces
      CGAL_assertion(is_insert_successful.first->second.size() == 1);
      is_insert_successful.first->second.push_back(icf);
    }
  }

  void convert_to_p2t2()
  {
    is_1_cover = true;

    p2t2.clear();
    p2t2.tds().set_dimension(2);

    cpp11::unordered_map<Vertex_handle /*dt2*/, Vertex_handle /*p2t2*/> vertex_correspondence_map;

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

    typedef boost::unordered_map<
              std::set<Vertex_handle>,
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
      fh->set_offsets(fit->vertex(0)->offset(),
                      fit->vertex(1)->offset(),
                      fit->vertex(2)->offset());

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

    std::cout << dt2.number_of_vertices() / 9 << " canonical vertices in dt2" << std::endl;
    std::cout << p2t2.number_of_vertices() << " vertices in p2t2" << std::endl;

    std::cout << cfc << " canonical faces in dt2" << std::endl;
    std::cout << p2t2.number_of_faces() << " canonical faces in p2t2" << std::endl;

    CGAL_postcondition(p2t2.tds().is_valid());
    CGAL_postcondition(p2t2.is_valid());
  }

  // number of canonical simplicies of each dimension
  size_type number_of_vertices() const
  {
    // @todo something better than this naive way (flag vertex/edge/face classes)
    size_type nv = 0;
    typename DT2::Finite_vertices_iterator vit = dt2.finite_vertices_begin(),
                                                                vend = dt2.finite_vertices_end();
    for(; vit!=vend; ++vit)
      if(is_canonical(vit))
        ++nv;

    return nv;
  }

//  size_type number_of_edges() const { }
//  size_type number_of_faces() const { }

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
  bool is_simplicial_complex() const { }

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

  bool is_canonical(const Face_handle fh) const
  {
    if(dt2.is_infinite(fh))
      return false;

    int min_off_x = std::numeric_limits<int>::max(),
        min_off_y = std::numeric_limits<int>::max();

    for(int i=0; i<3; i++)
    {
      Vertex_handle vh = fh->vertex(i);
      min_off_x = (std::min)(min_off_x, vh->offset().x());
      min_off_y = (std::min)(min_off_y, vh->offset().y());
    }

    return (min_off_x == 0 && min_off_y == 0);
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
    typename DT2::Faces_iterator fit = dt2.faces_begin(),
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

    typename DT2::Face_circulator fc = dt2.incident_faces(vh),
                                                       done = fc;
    do
    {
      fc->set_canonical_flag(is_canonical(fc));
    }
    while(++fc != done);
  }

  /// Low level functions to mascarade the DT2 as a periodic triangulation
  Offset compute_offset_shift(const Offset& off_1, const Offset& off_2) const
  {
    Offset shift_off((std::min)(off_1.x(), off_2.x()),
                     (std::min)(off_1.y(), off_2.y()));

    return shift_off;
  }

  // @todo something smarter, this function is core to everything else
  // How could we keep "canonical neighbors" in memory instead of having
  // to find them later...?
  Face_handle neighbor(Face_handle fh, int i) const
  {
    if(is_canonical(fh->neighbor(i)))
      return fh->neighbor(i);

    Vertex_handle e_vh_0 = fh->vertex((i+1)%3);
    Vertex_handle e_vh_1 = fh->vertex((i+2)%3);

    Offset shift_off = compute_offset_shift(e_vh_0->offset(), e_vh_1->offset());
    Offset ce_vh_0_off = e_vh_0->offset() - shift_off;

    Vertex_handle ce_vh_0 = periodic_vertices.at(
                              canonical_vertices.at(e_vh_0)).at(ce_vh_0_off);

    Vertex_handle cvh_1 = canonical_vertices.at(e_vh_1);

    Face_handle adj_fh;

    Face_circulator fc = dt2.incident_faces(ce_vh_0), done = fc;
    do
    {
      if(dt2.is_infinite(fc))
        continue;

      Vertex_handle ccw_ce_vh_0 = fc->vertex(dt2.ccw(fc->index(ce_vh_0)));

      if(canonical_vertices.at(ccw_ce_vh_0) != cvh_1)
        continue;

      if(ccw_ce_vh_0->offset() == (e_vh_1->offset() - shift_off))
        adj_fh = fc;
    }
    while(++fc != done);

    CGAL_assertion(adj_fh != Face_handle());

    // Now, we have the correct face, but it might not be canonical
    // @todo to find the canonical face, you can find the canonical offset
    // of the face, and then look at the difference of offset between the offset
    // of the vertices in the canonical face and in 'adj_fh'

    // @tmp
    Face_handle c_adj_fh = adj_fh;
    return c_adj_fh;

    CGAL_postcondition(false);
    return Face_handle();
  }

  /// Iterators and Circulators
  // @todo this should return a (custom) circulator
  std::set<Face_handle> incident_faces(Vertex_handle vh) const
  {
    std::set<Face_handle> ifhs;

    // Need to start from the canonical vertex to make sure there is no "weird" (to be defined properly) faces
    if(!is_canonical(vh))
      vh = canonical_vertices.at(vh);

    Face_circulator tds_fc = dt2.incident_faces(vh);

    // @todo check (when this function is switched to a Circulator if we must
    // walk cw or ccw)
    Face_handle fh = tds_fc->neighbor(dt2.ccw(tds_fc->index(vh))), done = fh;

    do
    {
      ifhs.insert(fh);

      // @todo when proper periodic traits are used to construct dt2 (that is,
      // insert (p, off) instead of constructing the offsetted point,
      // then the vertices will store the _canonical_ (geometric) point,
      // and we don't have to use the canonical_vertices[] but rather we can
      // just compare fh->vertex(0)->point() and vh->point()
      //
      // "Real" offsetted point is then done via a 'construct_point(cp, off)' function
      // in the triangulation class
      int vh_idx = -1;

      for(int i=0; i<3; ++i)
      {
        Vertex_handle vhi = fh->vertex(i);
        Vertex_handle cvhi = canonical_vertex(vhi);

        if(dt2.is_infinite(vhi))
          continue;

        // @fixme this is actually not sufficient, you might have the same vertex
        // appearing multiple times in the face. You need to look at the difference
        // of offsets between 'fh' and 'fhn' and deduce which vertex is the correct one.
        if(cvhi == vh)
          vh_idx = i;
      }

      CGAL_assertion(vh_idx != -1);

      fh = neighbor(fh, dt2.ccw(vh_idx));
    }
    while(fh != done);

    std::cout << ifhs.size() << " incident faces" << std::endl;

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
//  P2T2_Vertex_handle insert_into_p2t2(const Point& p)
//  {
//    const Point cp = lattice_.construct_canonical_point(p);
//    p2t2.insert(cp);
//  }

  Vertex_handle insert(const Point& p)
  {
    const Point cp = lattice_.construct_canonical_point(p);

    Vertex_handle vh = dt2.insert(cp);
    CGAL_assertion(vh != Vertex_handle());
    vh->set_offset(Offset(0, 0));

    mark_canonical_faces(vh);

    for(int off_x=-3; off_x<4; ++off_x)
    {
      for(int off_y=-3; off_y<4; ++off_y)
      {
        if(off_x == 0 && off_y == 0)
          continue;

        // @fixme Insert without constructions (cf P3T3)
        const Vector off_v = gt_.construct_sum_of_vectors_2_object()(
                               gt_.construct_scaled_vector_2_object()(
                                 lattice_.basis()[0], off_x),
                               gt_.construct_scaled_vector_2_object()(
                                 lattice_.basis()[1], off_y));

        const Point off_p = cp + off_v;

        if(lattice_.is_in_scaled_domain(off_p, 3))
        {
          Vertex_handle vh_copy = dt2.insert(off_p);
          CGAL_assertion(vh_copy != Vertex_handle());
          vh_copy->set_offset(Offset(off_x, off_y));

          canonical_vertices[vh_copy] = vh;
          periodic_vertices[vh][vh_copy->offset()] = vh_copy;

          mark_canonical_faces(vh_copy);
        }
      }
    }

    return vh;
  }

  // @todo should take two ForwardIterators
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

  DT2 dt2; // @tmp, shouldn't be exposed
  P2T2 p2t2;

private:
  Lattice lattice_; // @tmp to be incorporated in the geometric traits

  bool is_1_cover;

  cpp11::unordered_map<Vertex_handle /*periodic copy*/,
                       Vertex_handle /*canonical*/> canonical_vertices;

  // @todo hopefully don't need to use this eventually
  cpp11::unordered_map<Vertex_handle /*canonical*/,
                       std::map<Offset, // @todo unordered
                                Vertex_handle /*periodic copy*/> > periodic_vertices;

  Geom_traits gt_;
  Triangulation_data_structure tds_;
};

} //namespace CGAL

#endif // CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_GENERIC_H
