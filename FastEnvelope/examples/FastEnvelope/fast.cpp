
// #define CGAL_PROFILE


// #define TRACE

// fast.cpp and fastE.cpp produce the same trace output to see where the executables take different paths
// As fastE operates on reordered faces it is important to use as input an off file
// where the faces are ordered in such a way that they do not get reordered.
// Therefore fastE.cpp writess the indices into a file named "sorted.off" which can be used
// to replace the faces in the original input.


#include <Eigen/Dense>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_primitive.h>
#include <boost/iterator/counting_iterator.hpp>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <string>
#include <fstream>



namespace CGAL {
template <typename K>
int obtuse_angle(const Point_3<K>& p, const Point_3<K>& q, const Point_3<K>& r)
{
  if(angle(r,p,q) == OBTUSE){
    return 0;
  }
  if(angle(p,q,r) == OBTUSE){
    return 1;
  }
  if(angle(q,r,p) == OBTUSE){
    return 2;
  }
  return -1;
}

template <typename K>
Vector_3<K> normalize(const Vector_3<K>& v)
{
  return v / approximate_sqrt(v*v);
}


}// namespace CGAL


typedef Eigen::Matrix<int, 3, 1> Vector3i;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;



template <typename K>
struct Envelope {

  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Segment_3 Segment_3;
  typedef typename K::Triangle_3 Triangle_3;
  typedef typename K::Plane_3 Plane_3;
  typedef typename K::Iso_cuboid_3 Iso_cuboid_3;

  typedef CGAL::Bbox_3 Bbox_3;

  typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
  typedef typename EK::Point_3 ePoint_3;
  typedef typename EK::Line_3 eLine_3;
  typedef typename EK::Plane_3 ePlane_3;
  typedef typename EK::Intersect_3 eIntersect_3;
  typedef typename EK::Oriented_side_3 eOriented_side_3;
  typedef typename EK::Are_parallel_3 eAre_parallel_3;

  // The class  `Plane` is used for the 7-8 walls of a prism.
  // We store at the same  time threee points and a plane.
  // That is easier than retrieving the 3 points of a lazy plane.
  struct Plane {
    Plane()
    {}

    Plane(const ePoint_3& ep, const ePoint_3& eq, const ePoint_3& er)
      : ep(ep), eq(eq), er(er), eplane(ep,eq,er)
    {}

    Plane(const Point_3& p, const Point_3& q, const Point_3& r)
      : ep(p.x(),p.y(),p.z()), eq(q.x(),q.y(),q.z()), er(r.x(),r.y(),r.z()), eplane(ep,eq,er)
    {}
    ePoint_3 ep, eq, er;
    ePlane_3 eplane;
  };

  typedef std::vector<Plane> Prism;

  std::vector<Prism> halfspace; // should be renamed to "prisms"
  std::vector<Iso_cuboid_3> bounding_boxes;

  static const bool OUT_PRISM = 1;
  static const bool IN_PRISM = 0;
	static const int CUT_COPLANAR = 4;
	static const int CUT_EMPTY = -1;
	static const int CUT_FACE = 3;


  std::vector<Point_3> env_vertices;
  std::vector<Vector3i> env_faces;


  // For a query object the envelope test uses an AABB tree to find the relevant prisms

  // property maps for the primitive
  template <class Kernel>
  struct Datum_map
  {
    typedef boost::readable_property_map_tag category;
    typedef std::size_t key_type;
    typedef typename Kernel::Iso_cuboid_3 value_type;
    typedef value_type reference;

    const std::vector<Iso_cuboid_3>* boxes_ptr;

    Datum_map() : boxes_ptr(nullptr) {}
    Datum_map(const std::vector<Iso_cuboid_3>& boxes) : boxes_ptr(&boxes) {}

    friend value_type get(const Datum_map& m, key_type k)
    {
      return (*m.boxes_ptr)[k];
    }
  };

  template <class Kernel>
  struct Point_map
  {
    typedef boost::readable_property_map_tag category;
    typedef std::size_t key_type;
    typedef typename Kernel::Point_3 value_type;
    typedef value_type reference;

    const std::vector<Iso_cuboid_3>* boxes_ptr;

    Point_map() : boxes_ptr(nullptr) {}
    Point_map(const std::vector<Iso_cuboid_3>& boxes) : boxes_ptr(&boxes) {}

    friend value_type get(const Point_map& m, key_type k)
    {
      return ((*m.boxes_ptr)[k].min)();
    }
  };

  typedef CGAL::AABB_primitive<std::size_t, Datum_map<K>, Point_map<K>, CGAL::Tag_true /*UseSharedData*/, CGAL::Tag_false /*CacheDatum*/> Primitive;
  typedef CGAL::AABB_traits<K, Primitive> AABB_traits;
  typedef CGAL::AABB_tree<AABB_traits> Tree;


  Tree tree;

  eOriented_side_3 oriented_side;


  Envelope(const std::vector<Point_3>& env_vertices,
           std::vector<Vector3i> env_faces,
           double epsilon)
    : env_vertices(env_vertices), env_faces(env_faces)
  {
    halfspace_generation(env_vertices, env_faces, halfspace, bounding_boxes, epsilon);

    Datum_map<K> datum_map(bounding_boxes);
    Point_map<K> point_map(bounding_boxes);

    // constructs AABB tree
    tree.insert(boost::counting_iterator<std::size_t>(0),
                boost::counting_iterator<std::size_t>(bounding_boxes.size()),
                datum_map,
                point_map);
    tree.build();
  }


  struct INDEX
  {
    int Pi;
    std::vector<int> FACES;
  };


  // was marked todo???
  bool box_box_intersection(const Iso_cuboid_3 ic1, const Iso_cuboid_3 ic2) const
  {
    const Point_3& min1 = min_vertex(ic1);
    const Point_3& min2 = min_vertex(ic2);
    const Point_3& max1 = max_vertex(ic1);
    const Point_3& max2 = max_vertex(ic2);

    if ((compare_x(max1, min2) == CGAL::SMALLER) ||(compare_y(max1, min2) == CGAL::SMALLER) ||(compare_z(max1, min2) == CGAL::SMALLER)){
      return 0;
    }
    if ((compare_x(max2, min1) == CGAL::SMALLER) ||(compare_y(max2, min1) == CGAL::SMALLER) ||(compare_z(max2, min1) == CGAL::SMALLER)){
      return 0;
    }
    return 1;
  }


  bool do_intersect(const ePlane_3& ep0, const ePlane_3& ep1, const ePlane_3& ep2) const
  {

  CGAL::cpp11::result_of<eIntersect_3(ePlane_3, ePlane_3, ePlane_3)>::type
    result = CGAL::intersection(ep0,ep1,ep2);
  if(! result){
    return false;
  }

  const ePoint_3* ipp = boost::get<ePoint_3>(&*result);
  return ipp != nullptr;
  }


  bool
  point_out_prism(const ePoint_3 &point, const std::vector<unsigned int> &prismindex, int jump) const
  {
    //CGAL::Orientation ori;
    CGAL::Oriented_side ori;

    for (int i = 0; i < prismindex.size(); i++){
      if (prismindex[i] == jump){
        continue;
      }

      for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
        const Plane& plane = halfspace[prismindex[i]][j];
        ori = oriented_side(plane.eplane, point);
        //ori = CGAL::orientation(plane.ep, plane.eq, plane.er, point);
        if (ori != CGAL::ON_NEGATIVE_SIDE){
          // if for a prism we are on the wrong side of one halfspace we are outside this prism
          // so no need to look at the other halfspaces
          break;
        }
        if (j == halfspace[prismindex[i]].size() - 1){
          // As we are in all halfspaces of one prism we are in the union of the prisms
          return false;
        }
      }
    }

    return true;
  }


  // \param jump is the index of the prism that shall be ignored
  // \param id is a return parameter for the prism with `point` inside
  bool
  point_out_prism_return_local_id(const Point_3 &point, const ePoint_3 &epoint, const std::vector<unsigned int> &prismindex, const int jump, int &id) const
  {
    Vector_3 bmin, bmax;

    CGAL::Oriented_side ori;

    for (int i = 0; i < prismindex.size(); i++){
      if (prismindex[i] == jump){
        continue;
      }
      if(bounding_boxes[prismindex[i]].has_on_unbounded_side(point)){
        continue;
      }
      for (int j = 0; j < halfspace[prismindex[i]].size(); j++){
        const Plane& plane = halfspace[prismindex[i]][j];
        ori = oriented_side(plane.eplane, epoint);
        if (ori != CGAL::ON_NEGATIVE_SIDE){
          break;
        }

      }
      if (ori == CGAL::ON_NEGATIVE_SIDE){
        id = i;
        return false;
      }

    }

    return true;
  }



  // \param cindex  the index of a prism
  // \param cid     a return parameter where the indices of the faces that intersect the segment `(source,target)`get inserted
  bool
  is_seg_cut_polyhedra(const int cindex,
                       const ePoint_3& source,
                       const ePoint_3& target,
                       const eLine_3& line,
                       std::vector<int>& cid) const
  {
    const Prism& prism = halfspace[cindex];
    cid.clear();
    std::vector<bool> cut;

    cut.resize(prism.size());
    for (int i = 0; i < prism.size(); i++){
      cut[i] = false;
    }
    std::vector<CGAL::Oriented_side> o1, o2;
    o1.resize(prism.size());
    o2.resize(prism.size());
    int ori = 0, ct1 = 0, ct2 = 0;//ori=0 to avoid the case that there is only one cut plane
    std::vector<int> cutp;

    for (int i = 0; i < prism.size(); i++){
      const Plane& plane = prism[i];
      // POSITIVE is outside the prism
      o1[i] = oriented_side(plane.eplane, source);// CGAL::orientation(plane.ep, plane.eq, plane.er, source); // todo use plane.eplane
      o2[i] = oriented_side(plane.eplane, target);// CGAL::orientation(plane.ep, plane.eq, plane.er, target);

      if (int(o1[i]) + int(o2[i]) >= 1)
        {
          return false;
        }

      if (o1[i] == CGAL::ON_ORIENTED_BOUNDARY && o2[i] == CGAL::ON_ORIENTED_BOUNDARY)
        {
          return false;
        }

      if (int(o1[i]) * int(o2[i]) == -1){
        cutp.emplace_back(i);
      }
      if (o1[i] == CGAL::ON_POSITIVE_SIDE) ct1++;
      if (o2[i] == CGAL::ON_POSITIVE_SIDE) ct2++;// if ct1 or ct2 >0, then NOT totally inside
    }

    if (cutp.size() == 0 && ct1 == 0 && ct2 == 0){
      // no intersected planes, and each point is either inside of poly,
      //or on one facet, since vertices are checked, then totally inside
      return true;
    }
    if (cutp.size() == 0) {
      return false;
    }

    /* The segment can have an intersection with several planes,
       but they may be outside the prism.
       So we have to test for the intersection points i if they are outside the prism.



                                    |
                                    |
                                    |      t
                                    |     /
       -------------*****************----i-------------
                    *               *   /
                    *               *  /
                    *    prism      * /
                    *               ./
                    *               i
                    *              /.
                    *             / *
                    *            /  *
                    *           s   *
                    *               *
       -------------*****************-----------------
                    |               |
                    |               |
                    |               |
    */


    for (int i = 0; i < cutp.size(); i++){
      const Plane& plane_i = prism[cutp[i]];

      CGAL::cpp11::result_of<eIntersect_3(eLine_3, ePlane_3)>::type
        result = CGAL::intersection(line, plane_i.eplane);
      if(! result){
        std::cout <<  "there must be an intersection 2" << std::endl;
      }

      const ePoint_3* ipp = boost::get<ePoint_3>(&*result);
      CGAL_assertion(ipp != nullptr);

      const ePoint_3& ip = *ipp;

      for(int j = 0; j < cutp.size(); j++) {
        if (i == j){
          continue;
        }
        const Plane& plane_j = prism[cutp[j]];
        ori = oriented_side(plane_j.eplane, ip);

        if(ori == CGAL::ON_POSITIVE_SIDE){
          break;
        }
      }
      if (ori != CGAL::ON_POSITIVE_SIDE) {
        cut[cutp[i]] = true;
      }
    }

    for (int i = 0; i < prism.size(); i++) {
      if (cut[i] == true){
        cid.emplace_back(i);
      }
    }
    if (cid.size() == 0){
      return false;// if no intersection points, and segment not totally inside, then not intersected
    }
    return true;
  }


  int
  Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(const ePoint_3 &source, const ePoint_3 &target,
                                                          const eLine_3& eline,
                                                          const Plane& plane,
                                                          const std::vector<unsigned int> &prismindex, const int &jump, int &id) const
  {
    CGAL::Oriented_side ori;
    CGAL::cpp11::result_of<eIntersect_3(eLine_3, ePlane_3)>::type
      result = CGAL::intersection(eline, plane.eplane);
#ifdef TRACE
    if(! result){
      std::cout <<  "there must be an intersection 3" << std::endl;
    }
#endif

    const ePoint_3* ipp = boost::get<ePoint_3>(&*result);
    CGAL_assertion(ipp != nullptr);

    const ePoint_3& ip = *ipp;
    for (int i = 0; i < prismindex.size(); i++){
      if (prismindex[i] == jump){
        continue;
      }

      for (int j = 0; j < halfspace[prismindex[i]].size(); j++){
        const Plane& plane = halfspace[prismindex[i]][j];
        ori = oriented_side(plane.eplane, ip);

        if (ori != CGAL::ON_NEGATIVE_SIDE){
          break;
        }
      }
      if (ori == CGAL::ON_NEGATIVE_SIDE){
        id = i;
        return IN_PRISM;
      }
    }

    return OUT_PRISM;

  }



  bool
  segment_out_of_envelope(const Point_3& source, const Point_3& target,
                          const std::vector<unsigned int>& prismindex) const
  {
    if (prismindex.size() == 0) {
      return true;
    }


    ePoint_3 esource(source.x(),source.y(),source.z());
    ePoint_3 etarget(target.x(),target.y(),target.z());

    int jump1 = -1;
    bool out, cut;
    int inter;

    int check_id, check_id1;

    // First check if the endpoints are outside the envelope
    out = point_out_prism_return_local_id(source, esource, prismindex, jump1, check_id);

    if (out) {
      return true;
    }
    out = point_out_prism_return_local_id(target, etarget, prismindex, jump1, check_id1);

    if (out) {
      return true;
    }

    // If both endpoints are in the same prism it is in the envelope
    if (check_id == check_id1){
      return false;
    }
    if (prismindex.size() == 1){
      return false;
    }
    eLine_3 line(esource,etarget);
    std::vector<unsigned int > queue, idlist;
    queue.emplace_back(check_id);//queue contains the id in prismindex
    idlist.emplace_back(prismindex[check_id]);

    std::vector<int> cidl;
    cidl.reserve(8);
    for (int i = 0; i < queue.size(); i++) {

      jump1 = prismindex[queue[i]];

      cut = is_seg_cut_polyhedra(jump1, esource, etarget, line, cidl);
      // cidl now contains the faces of the prism jump1
      if (cut&&cidl.size() == 0){
        return false;
      }
      if (!cut){
        continue;
      }

      for (int j = 0; j < cidl.size(); j++) {

        inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id
          (esource,
           etarget,
           line,
           halfspace[prismindex[queue[i]]][cidl[j]],
           idlist, jump1, check_id);

        if (inter == 1){
          inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id
            (esource,
             etarget,
             line,
             halfspace[prismindex[queue[i]]][cidl[j]],
             prismindex, jump1, check_id);

          if (inter == 1) {
            return true; // outside envelope
          }
          if (inter == 0) {
            queue.emplace_back(check_id);
            idlist.emplace_back(prismindex[check_id]);
          }
        }
      }
    }

    return false; // fully inside the envelope
  }



  int
  is_3_triangle_cut_float_fast(const ePoint_3& tp,
                               const ePoint_3& tq,
                               const ePoint_3& tr,
                               const ePlane_3 &tri,
                               const ePlane_3 &facet1,
                               const ePlane_3 &facet2) const
  {
    // todo:  what do we test here with n ?
    // todo : do this not with Epeck
    ePoint_3 n = tp + CGAL::cross_product((tp - tq), (tp - tr));

    if (CGAL::orientation(n, tp, tq, tr) == 0){
        std::cout << "todo degeneration handling" << std::endl;
        //n = Point_3(rand(), rand(), rand())} };
      }

    CGAL::cpp11::result_of<eIntersect_3(ePlane_3, ePlane_3, ePlane_3)>::type
      result = CGAL::intersection(tri, facet1, facet2);
    if(! result){
#ifdef TRACE
      std::cout <<  "there must be an intersection 4" << std::endl;
#endif
        return 0;
    }
    const ePoint_3* ipp = boost::get<ePoint_3>(&*result);
    CGAL_assertion(ipp != nullptr);

    const ePoint_3& ip = *ipp;

    int o1 = int(CGAL::orientation(n,tp,tq, ip));
    int o2 = int(CGAL::orientation(n,tq,tr, ip));

    if (o1 * o2 == -1){
      return 0;
    }
    int o3 = int(CGAL::orientation(n, tr, tp, ip));

    if (o1 * o3 == -1 || o2 * o3 == -1)
      return 0;
    if (o1 * o2 * o3 == 0)
      return 2; // means we dont know
    return 1;
  }


  bool
  is_3_triangle_cut(const ePoint_3& tp,
                    const ePoint_3& tq,
                    const ePoint_3& tr,
                    const ePlane_3 &tri,
                    const ePlane_3 &facet1,
                    const ePlane_3 &facet2) const
  {

    CGAL::cpp11::result_of<eIntersect_3(ePlane_3, ePlane_3, ePlane_3)>::type
      result = CGAL::intersection(tri, facet1, facet2);
    if(! result){
      // todo:  what to do?
#ifdef TRACE
      std::cout <<  "there must be an intersection 5" << std::endl;
#endif
    }
    const ePoint_3* ipp = boost::get<ePoint_3>(&*result);
    if(ipp != nullptr){
      // todo:  what to do?
#ifdef TRACE
      std::cout <<  "the intersection must be a point" << std::endl;
#endif
    }

    const ePoint_3& ip = *ipp;


    ePoint_3 n = tp + CGAL::cross_product((tp - tq), (tp - tr));

    if (Predicates::orient_3d(n, tp, q, tr) == 0)
      {
        //logger().debug("Degeneration happens");
#ifdef TRACE
        std::cout << "todo degeneration handling" << std::endl;
#endif
        //n = { {Vector3(rand(), rand(), rand())} };
      }

    int o1 = int(orient(n,tp,tq, ip));
    if (o1 == 0){
      return false;
    }

    int o2 = int(orientation(n, tq, tr, ip));

    if (o2 == 0 || o1 + o2 == 0){
      return false;
    }

    int o3 = int(orientation(n, tr, tp, ip));

    if (o3 == 0 || o1 + o3 == 0 || o2 + o3 == 0){
      return false;
    }

    return true;
  }


  bool
  is_two_facets_neighbouring(const int & pid, const int &i, const int &j)const
  {
    int facesize = halfspace[pid].size();
    if (i == j) return false;
    if (i == 0 && j != 1) return true;
    if (i == 1 && j != 0) return true;
    if (j == 0 && i != 1) return true;
    if (j == 1 && i != 0) return true;
    if (i - j == 1 || j - i == 1) return true;
    if (i == 2 && j == facesize - 1) return true;
    if (j == 2 && i == facesize - 1) return true;
    return false;
  }


  int
  is_triangle_cut_envelope_polyhedra(const int &cindex,//the triangle is not degenerated
                                     const ePoint_3 &tri0, const ePoint_3 &tri1, const ePoint_3 &tri2,
                                     std::vector<int> &cid) const
  {
    const Prism& prism = halfspace[cindex];
    cid.clear();
    cid.reserve(3);
    std::vector<bool> cut;
    cut.resize(prism.size());
    for (int i = 0; i < prism.size(); i++)
      {
        cut[i] = false;
      }
    std::vector<int> o1, o2, o3, cutp;
    o1.resize(prism.size());
    o2.resize(prism.size());
    o3.resize(prism.size());
    int  ori = 0, ct1 = 0, ct2 = 0, ct3 = 0;


    for (int i = 0; i < prism.size(); i++)
      {
        const Plane& plane = prism[i];
        o1[i] = int(oriented_side(plane.eplane, tri0));
        o2[i] = int(oriented_side(plane.eplane, tri1));
        o3[i] = int(oriented_side(plane.eplane, tri2));
        if (o1[i] + o2[i] + o3[i] >= 2){ //1,1,0 case
          return 0;
        }
        if (o1[i] == 1) ct1++;
        if (o2[i] == 1) ct2++;
        if (o3[i] == 1) ct3++;// if ct1 or ct2 or ct3 >0, then NOT totally inside, otherwise, totally inside
        if (o1[i] == 0 && o2[i] == 0 && o3[i] == 1){
            return 0;
          }
        if (o1[i] == 1 && o2[i] == 0 && o3[i] == 0){
            return 0;
          }
        if (o1[i] == 0 && o2[i] == 1 && o3[i] == 0){
            return 0;
          }


        if (o1[i] * o2[i] == -1 || o1[i] * o3[i] == -1 || o3[i] * o2[i] == -1){
          cutp.emplace_back(i);
        }else if (o1[i] + o2[i] + o3[i] == -1 && o1[i] * o2[i] == 0) {//0,0,-1 case, we also want this face,really rare to happen
          cutp.emplace_back(i);
        }
      }
    if (cutp.size() == 0) {
      if (ct1 == 0 && ct2 == 0 && ct3 == 0) {
        return 2;// totally inside, or not any edge is on the facet
      }
    }

    if (cutp.size() == 0){
        return 0;
      }
#ifdef TRACE
    std::cout << "A" << std::endl;
#endif
    std::array<ePoint_3*,2> seg0, seg1;

    for (int i = 0; i < cutp.size(); i++)
      {
        int tmp = 0;
        if (o1[cutp[i]] * o2[cutp[i]] == -1|| o1[cutp[i]] + o2[cutp[i]] == -1) {
          seg0[tmp] = const_cast<ePoint_3*>(&tri0);
          seg1[tmp] = const_cast<ePoint_3*>(&tri1);

          tmp++;
        }
        if (o1[cutp[i]] * o3[cutp[i]] == -1|| o1[cutp[i]] + o3[cutp[i]] == -1) {
          seg0[tmp] = const_cast<ePoint_3*>(&tri0);
          seg1[tmp] = const_cast<ePoint_3*>(&tri2);

          tmp++;
        }
        if (o2[cutp[i]] * o3[cutp[i]] == -1|| o2[cutp[i]] + o3[cutp[i]] == -1) {
          seg0[tmp] = const_cast<ePoint_3*>(&tri1);
          seg1[tmp] = const_cast<ePoint_3*>(&tri2);

          tmp++;
        }


        for (int k = 0; k < 2; k++){
          const Plane& plane_i = prism[cutp[i]];

          //         std::cout <<  plane_i.ep << "  " <<  plane_i.eq << "  " <<  plane_i.er << "  " <<  plane_i.ep << std::endl;
          // std::cout << *(seg0[k]) << "       " <<   *(seg1[k]) << std::endl;

          eLine_3 eline(*(seg0[k]), *(seg1[k]));
          CGAL::cpp11::result_of<eIntersect_3(eLine_3, ePlane_3)>::type
            result = CGAL::intersection(eline, plane_i.eplane);
          if(! result){
#ifdef TRACE
            std::cout <<  "there must be an intersection 6" << std::endl;
#endif
          }

          const ePoint_3* ipp = boost::get<ePoint_3>(&*result);
          CGAL_assertion(ipp != nullptr);
          const ePoint_3& ip = *ipp;

          for (int j = 0; j < cutp.size(); j++){
            if (i == j){
                  continue;
            }
            const Plane& plane_j = prism[cutp[j]];
            ori = oriented_side(plane_j.eplane, ip); //  todo:  check  if it must be positive or negative
            if (ori == 1){
              break;
            }
          }
          if (ori != 1){
            int cutpi = cutp[i];
            cut[cutp[i]] = true;
            break;
          }
        }

        ori = 0;// initialize the orientation to avoid the j loop doesn't happen because cutp.size()==1
      }

#ifdef TRACE
    std::cout << "B" << std::endl;
#endif
    if (cutp.size() <= 2){
        for (int i = 0; i < prism.size(); i++){
          if (cut[i] == true){
              cid.emplace_back(i);
          }
        }
        return 1;
    }

#ifdef TRACE
      std::cout << "C" << std::endl;
#endif
    // triangle-facet-facet intersection

    ePlane_3 tri_eplane(tri0, tri1, tri2); // todo change to a query triangle with 3 points and a plane
    for (int i = 0; i < cutp.size(); i++)
      {
        for (int j = i + 1; j < cutp.size(); j++)// two facets and the triangle generate a point
          {
#ifdef TRACE
            std::cout << "T :\n "<< tri0 << "    "  << tri1 << "    " << tri2 << std::endl;
            std::cout << "I : "<< i << " " << cutp[i] << "\n" <<  prism[cutp[i]].ep << std::endl <<  prism[cutp[i]].eq << std::endl <<  prism[cutp[i]].er << std::endl;
            std::cout << "J : "<< j << " " << cutp[j] << "\n" <<  prism[cutp[j]].ep << std::endl <<  prism[cutp[j]].eq << std::endl <<  prism[cutp[j]].er << std::endl;
#endif
            if (cut[cutp[i]] == true && cut[cutp[j]] == true)
              continue;

            if (true /* USE_ADJACENT_INFORMATION*/ ) {
              bool neib = is_two_facets_neighbouring(cindex, cutp[i], cutp[j]);
              if (neib == false) continue;
            }


            //            bool inter = this->do_intersect(tri_eplane,prism[cutp[i]].eplane, prism[cutp[j]].eplane);
            int inter = is_3_triangle_cut_float_fast(tri0, tri1,tri2,tri_eplane,prism[cutp[i]].eplane, prism[cutp[j]].eplane);
#ifdef TRACE
            std::cout << "is_3_triangle_cut_float_fast: " << inter << std::endl;
#endif
            // this was for a fast float check
            if (inter == 2)
              { //we dont know if point exist or if inside of triangle
                int cutpi = cutp[i];
                int cutpj = cutp[j];
                cut[cutp[i]] = true;
                cut[cutp[j]] = true;
                continue;
              }

            if (inter == 0){
              continue; // sure not inside
            }

            CGAL::cpp11::result_of<eIntersect_3(ePlane_3, ePlane_3, ePlane_3)>::type
              result = CGAL::intersection(tri_eplane, prism[cutp[i]].eplane, prism[cutp[j]].eplane);
            if(! result){
#ifdef TRACE
              std::cout <<  "there must be an intersection 7" << std::endl;
#endif
            }
            const ePoint_3* ipp = boost::get<ePoint_3>(&*result);
            CGAL_assertion(ipp != nullptr);

            const ePoint_3& ip = *ipp;


            for (int k = 0; k < cutp.size(); k++){

              if (k == i || k == j){
                  continue;
              }
#ifdef TRACE
              std::cout << k << " " << cutp[k] << "\n" <<  prism[cutp[k]].ep << std::endl <<  prism[cutp[k]].eq << std::endl <<  prism[cutp[k]].er << std::endl;
#endif
              ori = int(oriented_side(prism[cutp[k]].eplane, ip));
#ifdef TRACE
              std::cout << ori << std::endl;
#endif

              if (ori == 1){
                  break;
              }
            }

            if (ori != 1) {
              int cutpi = cutp[i];
              int cutpj = cutp[j];
              cut[cutp[i]] = true;
              cut[cutp[j]] = true;
            }
          }
      }
#ifdef TRACE
    std::cout << "D" << std::endl;
#endif

    for (int i = 0; i < prism.size(); i++){
      if (cut[i] == true){
#ifdef TRACE
             std::cout << "cut " << i << " is true" << std::endl;
#endif
          cid.emplace_back(i);
      }
      }
    //    std::cout << "cid.size()= " << cid.size() << std::endl;
    if (cid.size() == 0){
      return 0;// not cut and facets, and not totally inside, then not intersected
    }
    return 1;
  }


  int
  seg_cut_plane(const ePoint_3 &seg0, const ePoint_3 &seg1,
                const ePoint_3 &t0, const ePoint_3 &t1, const ePoint_3 &t2) const
  {
    int o1, o2;
    o1 = int(CGAL::orientation(seg0, t0, t1, t2));
    o2 = int(CGAL::orientation(seg1, t0, t1, t2));
    int op = o1 * o2;
    if (op >= 0)
      {
        return CUT_COPLANAR; //in fact, coplanar and not cut this plane
      }
    return CUT_FACE;
  }


  bool
  is_tpp_on_polyhedra(const ePlane_3 &triangle,
                      const ePlane_3 &facet1,
                      const ePlane_3 &facet2,
                      const int &prismid, const int &faceid)const
  {
      CGAL::cpp11::result_of<eIntersect_3(ePlane_3, ePlane_3, ePlane_3)>::type
      result = CGAL::intersection(triangle, facet1, facet2);
      if(! result){
#ifdef TRACE
        std::cout <<  "there must be an intersection 8" << std::endl;
#endif
      }

      const ePoint_3* ipp = boost::get<ePoint_3>(&*result);
      CGAL_assertion(ipp != nullptr);

      const ePoint_3& ip = *ipp;

       for (int i = 0; i < halfspace[prismid].size(); i++) {
        /*bool neib = is_two_facets_neighbouring(prismid, i, faceid);// this works only when the polyhedron is convex and no two neighbour facets are coplanar
          if (neib == false) continue;*/
        if (i == faceid) continue;
        if(oriented_side(halfspace[prismid][i].eplane, ip) == CGAL::ON_POSITIVE_SIDE){ //todo check positive/negative
          return false;
        }
       }
       return true;
  }



  int
  Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order(
		const ePoint_3 &segpoint0, const ePoint_3 &segpoint1,
                const ePlane_3 &eplane,
                const std::vector<unsigned int> &prismindex,
		const std::vector<std::vector<int>>& intersect_face, const int &jump, int &id) const
  {
    eLine_3 eline(segpoint0,segpoint1); // todo replace parameter of function

    CGAL::cpp11::result_of<eIntersect_3(eLine_3, ePlane_3)>::type
      result = CGAL::intersection(eline, eplane);
    if(! result){
 #ifdef TRACE
      std::cout <<  "there must be an intersection 9" << std::endl;
 #endif
    }

    const ePoint_3* ipp = boost::get<ePoint_3>(&*result);
    CGAL_assertion(ipp != nullptr);

    const ePoint_3& ip = *ipp;
    int tot, fid, ori;
    for (int i = 0; i < prismindex.size(); i++){
      if (prismindex[i] == jump){
          continue;
        }
      tot = 0; fid = 0;
      ori = -1;
      const Prism& prism = halfspace[prismindex[i]];
      for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {

        if (intersect_face[i][fid] == j)
          {
            if (fid + 1 < intersect_face[i].size()) fid++;
          }
        else continue;
        ori = int(oriented_side(prism[j].eplane, ip));

        if (ori == 1 || ori == 0)
          {
            break;
          }

        if (ori == -1)
          {
            tot++;
          }
      }
      if (ori == 1 || ori == 0) continue;
      fid = 0;
      ori = -1;
      for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
        if (intersect_face[i][fid] == j)
          {
            if (fid + 1 < intersect_face[i].size()) fid++;
            continue;
          }

        ori = int(oriented_side(prism[j].eplane, ip));
        if (ori == 1 || ori == 0)
          {
            break;
          }

        if (ori == -1)
          {
            tot++;
          }
      }
      if (ori == 1 || ori == 0) continue;
      if (tot == halfspace[prismindex[i]].size())
        {
          id = i;
          return IN_PRISM;
        }
    }


      return OUT_PRISM;
  }


  int
  Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over(
		const ePoint_3 &segpoint0, const ePoint_3 &segpoint1,
                const ePlane_3 &eplane,
                const std::vector<unsigned int> &prismindex,
		const std::vector<std::vector<int>>& intersect_face,
                const std::vector<bool>& coverlist,
                const int &jump,
                int &id) const
  {
    eLine_3 eline(segpoint0,segpoint1); // todo replace parameter of function

    CGAL::cpp11::result_of<eIntersect_3(eLine_3, ePlane_3)>::type
      result = CGAL::intersection(eline, eplane);
    if(! result){
#ifdef TRACE
      std::cout <<  "there must be an intersection 10" << std::endl;
#endif
    }

    const ePoint_3* ipp = boost::get<ePoint_3>(&*result);
    CGAL_assertion(ipp != nullptr);

    const ePoint_3& ip = *ipp;

    int tot, ori, fid;

    for (int i = 0; i < prismindex.size(); i++){
      if (prismindex[i] == jump){
        continue;
      }
      if (coverlist[i] == true){
        continue;
      }
      tot = 0; fid = 0;
      ori = -1;
      const Prism& prism = halfspace[prismindex[i]];

      for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
        if (intersect_face[i][fid] == j)   {
          if (fid + 1 < intersect_face[i].size()) fid++;
        }else{
          continue;
        }
        ori = int(oriented_side(prism[j].eplane,ip));
        if (ori == 1 || ori == 0){
          break;
        }

        if (ori == -1){
          tot++;
        }
      }
      if (ori == 1 || ori == 0){
        continue;
      }
      fid = 0;
      ori = -1;

      for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
        if (intersect_face[i][fid] == j){
          if (fid + 1 < intersect_face[i].size()) fid++;{
            continue;
          }
        }

        ori = int(oriented_side(prism[j].eplane,ip));
        if (ori == 1 || ori == 0){
          break;
        }

        if (ori == -1){
          tot++;
        }
      }
      if (ori == 1 || ori == 0){
        continue;
      }
      if (tot == halfspace[prismindex[i]].size()) {
        id = i;
        return IN_PRISM;
      }
    }

    return OUT_PRISM;
  }


  int
  Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order(
		const ePlane_3& triangle,
		const ePlane_3& facet1,
                const ePlane_3& facet2,
		const std::vector<unsigned int> &prismindex,
                const std::vector<std::vector<int>>&intersect_face,
                const int &jump1,
                const int &jump2,
		int &id) const
  {
      CGAL::cpp11::result_of<eIntersect_3(ePlane_3, ePlane_3, ePlane_3)>::type
      result = CGAL::intersection(triangle, facet1, facet2);
      if(! result){
 #ifdef TRACE
        std::cout <<  "there must be an intersection 11" << std::endl;
 #endif
      }

      const ePoint_3* ipp = boost::get<ePoint_3>(&*result);
      CGAL_assertion(ipp != nullptr);

      const ePoint_3& ip = *ipp;
      int tot, ori, fid;
      for (int i = 0; i < prismindex.size(); i++)
        {

          if (prismindex[i] == jump1 || prismindex[i] == jump2)	continue;
          if (!box_box_intersection(bounding_boxes[prismindex[i]], bounding_boxes[jump1])) continue;
          if (!box_box_intersection(bounding_boxes[prismindex[i]], bounding_boxes[jump2])) continue;

          tot = 0;
          fid = 0;
          ori = -1;
          const Prism& prism = halfspace[prismindex[i]];
          for (int j = 0; j < prism.size(); j++) {
            if (intersect_face[i][fid] == j)
              {
                if (fid + 1 < intersect_face[i].size()) fid++;
              }
            else continue;

            ori = int(oriented_side(prism[j].eplane, ip));

            if (ori == 1 || ori == 0)
              {
                break;
              }

            if (ori == -1)
              {
                tot++;
              }
          }
          if (ori == 1 || ori == 0) continue;
          fid = 0;
          ori = -1;
          for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            if (intersect_face[i][fid] == j)
              {
                if (fid + 1 < intersect_face[i].size()) fid++;
                continue;
              }


            ori = int(oriented_side(prism[j].eplane, ip));

            if (ori == 1 || ori == 0)
              {
                break;
              }

            if (ori == -1)
              {
                tot++;
              }
          }
          if (ori == 1 || ori == 0) continue;
          if (tot == prism.size())
            {
              id = i;
              return IN_PRISM;
            }

        }

      return OUT_PRISM;
    }



  int
  Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over(
		const ePlane_3 &triangle,
		const ePlane_3 &facet1,
                const ePlane_3 &facet2,
		const std::vector<unsigned int>& prismindex,
                const std::vector<std::vector<int>>& intersect_face,
                const std::vector<bool>& coverlist,
                const int &jump1,
                const int &jump2,
		int &id) const
  {
      CGAL::cpp11::result_of<eIntersect_3(ePlane_3, ePlane_3, ePlane_3)>::type
      result = CGAL::intersection(triangle, facet1, facet2);
      if(! result){
 #ifdef TRACE
        std::cout <<  "there must be an intersection 12" << std::endl;
#endif
      }

      const ePoint_3* ipp = boost::get<ePoint_3>(&*result);
      CGAL_assertion(ipp != nullptr);

      const ePoint_3& ip = *ipp;

      int tot, ori, fid;
      for (int i = 0; i < prismindex.size(); i++)
        {

          if (prismindex[i] == jump1 || prismindex[i] == jump2)	continue;
          if (!box_box_intersection(bounding_boxes[prismindex[i]], bounding_boxes[jump1])) continue;
          if (!box_box_intersection(bounding_boxes[prismindex[i]], bounding_boxes[jump2])) continue;
          if (coverlist[i] == true) continue;
          tot = 0;
          fid = 0;
          ori = -1;
          const Prism& prism = halfspace[prismindex[i]];
          for (int j = 0; j < prism.size(); j++) {
            if (intersect_face[i][fid] == j)
              {
                if (fid + 1 < intersect_face[i].size()) fid++;
              }
            else continue;

            ori = int(oriented_side(prism[j].eplane, ip));

            if (ori == 1 || ori == 0)
              {
                break;
              }

            if (ori == -1)
              {
                tot++;
              }
          }
          if (ori == 1 || ori == 0) continue;
          fid = 0;
          ori = -1;
          for (int j = 0; j < halfspace[prismindex[i]].size(); j++) {
            if (intersect_face[i][fid] == j)
              {
                if (fid + 1 < intersect_face[i].size()) fid++;
                continue;
              }

            ori = int(oriented_side(prism[j].eplane, ip));

            if (ori == 1 || ori == 0)
              {
                break;
              }

            if (ori == -1)
              {
                tot++;
              }
          }
          if (ori == 1 || ori == 0) continue;
          if (tot == halfspace[prismindex[i]].size())
            {
              id = i;
              return IN_PRISM;
            }

        }

      return OUT_PRISM;
    }



  bool
  triangle_out_of_envelope(const std::array<Point_3,3>& triangle,
                           const std::vector<unsigned int> &prismindex) const
  {
    if (prismindex.size() == 0)
      {
        return true;
      }

    std::array<ePoint_3,3> etriangle = { ePoint_3(triangle[0].x(),triangle[0].y(),triangle[0].z()),
                                         ePoint_3(triangle[1].x(),triangle[1].y(),triangle[1].z()),
                                         ePoint_3(triangle[2].x(),triangle[2].y(),triangle[2].z()) };

    int jump1, jump2;
    static const std::array<std::array<int, 2>, 3> triseg = {
      {{{0, 1}}, {{0, 2}}, {{1, 2}}}
    };


    std::vector<unsigned int> filtered_intersection;
    filtered_intersection.reserve(prismindex.size() / 3);
    std::vector<std::vector<int>>intersect_face;
    intersect_face.reserve(prismindex.size() / 3);
    bool out, cut;

    int inter, inter1, record1, record2,

      tti; //triangle-triangle intersection

    jump1 = -1;

    int check_id;

    for (int i = 0; i < 3; i++) {
      out = point_out_prism_return_local_id(triangle[i], etriangle[i], prismindex, jump1, check_id);

      if (out) {
        return true;
      }
    }

    if (prismindex.size() == 1)
      return false;


#ifdef DEGENERATION_FIX

    int degeneration = algorithms::is_triangle_degenerated(triangle[0], triangle[1], triangle[2]);

    if (degeneration == DEGENERATED_POINT){ //case 1 degenerate to a point
      return false;
    }

    if (degeneration == DEGENERATED_SEGMENT){
      std::vector<unsigned int > queue, idlist;
      queue.emplace_back(check_id);//queue contains the id in prismindex
      idlist.emplace_back(prismindex[check_id]);

      std::vector<int> cidl; cidl.reserve(8);
      for (int i = 0; i < queue.size(); i++) {

        jump1 = prismindex[queue[i]];
        int seg_inside = 0;
        for (int k = 0; k < 3; k++) {
          eLine_3 eline(etriangle[triseg[k][0]], etriangle[triseg[k][1]]);  // todo: store 3 lines in a triangle query object
          cut = is_seg_cut_polyhedra(jump1, etriangle[triseg[k][0]], etriangle[triseg[k][1]], eline, cidl);
          if (cut&&cidl.size() == 0) {
            seg_inside++;
            if (seg_inside == 3) return false;// 3 segs are all totally inside of some polyhedrons
            continue;// means this seg is inside, check next seg
          }
          if (!cut) continue;

          for (int j = 0; j < cidl.size(); j++) {
            inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(etriangle[triseg[k][0]], etriangle[triseg[k][1]],
                                                                            eline,
                                                                            halfspace[prismindex[queue[i]]][cidl[j]],
                                                                            idlist, jump1, check_id);


            if (inter == 1){
                inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id(etriangle[triseg[k][0]], etriangle[triseg[k][1]],
                                                                                eline,
                                                                                halfspace[prismindex[queue[i]]][cidl[j]],
                                                                                prismindex, jump1, check_id);



                if (inter == 1) {
                  return true;
                }
                if (inter == 0) {
                  queue.emplace_back(check_id);
                  idlist.emplace_back(prismindex[check_id]);
                }
              }
          }
        }
      }

      return false;
    }
    //
#endif // DEGENERATION_FIX


    std::vector<int> cidl; cidl.reserve(8);
    for (int i = 0; i < prismindex.size(); i++) {
      tti = is_triangle_cut_envelope_polyhedra(prismindex[i],
                                               etriangle[0], etriangle[1], etriangle[2], cidl);
      if (tti == 2) {

        return false;//totally inside of this polyhedron
      }
      else if (tti == 1 && cidl.size() > 0) {
        filtered_intersection.emplace_back(prismindex[i]);
        intersect_face.emplace_back(cidl);


      }
    }

    if (filtered_intersection.size() == 0) {
      return false;
    }

    //    std::cout << filtered_intersection.size() << " filtered" << std::endl;

    for(int i = 0; i < filtered_intersection.size(); i++){
      prism_to_off(filtered_intersection[i], "filtered");
    }

    std::vector<unsigned int > queue, idlist;
    std::vector<bool> coverlist;
    coverlist.resize(filtered_intersection.size());
    for (int i = 0; i < coverlist.size(); i++) {
      coverlist[i] = false;// coverlist shows if the element in filtered_intersection is one of the current covers
    }
    queue.emplace_back(0);//queue contains the id in filtered_intersection
    idlist.emplace_back(filtered_intersection[queue[0]]);// idlist contains the id in prismid//it is fine maybe it is not really intersected
    coverlist[queue[0]] = true;//when filtered_intersection[i] is already in the cover list, coverlist[i]=true

    std::vector<unsigned int> neighbours;//local id
    std::vector<unsigned int > list;
    std::vector<std::vector<int>> neighbour_facets, idlistorder;
    std::vector<bool> neighbour_cover;
    idlistorder.emplace_back(intersect_face[queue[0]]);

    for (int i = 0; i < queue.size(); i++) {

      jump1 = filtered_intersection[queue[i]];

      for (int k = 0; k < 3; k++) {
        for (int j = 0; j < intersect_face[queue[i]].size(); j++) {
          tti = seg_cut_plane(etriangle[triseg[k][0]],
                                          etriangle[triseg[k][1]],
                                          halfspace[filtered_intersection[queue[i]]][intersect_face[queue[i]][j]].ep,
                                          halfspace[filtered_intersection[queue[i]]][intersect_face[queue[i]][j]].eq,
                                          halfspace[filtered_intersection[queue[i]]][intersect_face[queue[i]][j]].er);

          if (tti != CUT_FACE) continue;

          inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order(etriangle[triseg[k][0]],
                                                                                          etriangle[triseg[k][1]],
                                                                                          halfspace[filtered_intersection[queue[i]]][intersect_face[queue[i]][j]].eplane,
                                                                                          idlist, idlistorder, jump1, check_id);



          if (inter == 1)
            {

              inter = Implicit_Seg_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over(etriangle[triseg[k][0]], etriangle[triseg[k][1]],
                                                                                                        halfspace[filtered_intersection[queue[i]]][intersect_face[queue[i]][j]].eplane,
                                                                                                        filtered_intersection, intersect_face, coverlist, jump1, check_id);


              assert(inter != 2);// the point must exist because it is a seg-halfplane intersection
              if (inter == 1) {

                return true;
              }
              if (inter == 0) {
                idlistorder.emplace_back(intersect_face[check_id]);
                queue.emplace_back(check_id);
                idlist.emplace_back(filtered_intersection[check_id]);
                coverlist[check_id] = true;
              }
            }
        }
      }
    }




    //tpi part

    //tree

    Tree localtree;

    std::vector<Iso_cuboid_3> local_bounding_boxes;
    local_bounding_boxes.resize(filtered_intersection.size());

    for (int i = 0; i < filtered_intersection.size(); i++) {
      local_bounding_boxes[i] = bounding_boxes[filtered_intersection[i]];
    }

    Datum_map<K> datum_map(local_bounding_boxes);
    Point_map<K> point_map(local_bounding_boxes);

    // constructs AABB tree
    localtree.insert(boost::counting_iterator<std::size_t>(0),
                      boost::counting_iterator<std::size_t>(local_bounding_boxes.size()),
                      datum_map,
                      point_map);
    localtree.build();

    //tree end

    for (int i = 1; i < queue.size(); i++){
      jump1 = filtered_intersection[queue[i]];

      localtree.all_intersected_primitives(bounding_boxes[jump1], std::back_inserter(list));
      neighbours.clear();
      neighbour_cover.clear();
      neighbour_facets.clear();

      neighbours.resize(list.size());
      neighbour_facets.resize(list.size());
      neighbour_cover.resize(list.size());
      for (int j = 0; j < list.size(); j++) {
        neighbours[j] = filtered_intersection[list[j]];
        neighbour_facets[j] = intersect_face[list[j]];
        if (coverlist[list[j]] == true) neighbour_cover[j] = true;
        else neighbour_cover[j] = false;
      }


      ePlane_3 etriangle_eplane(etriangle[0],etriangle[1],etriangle[2]);

      for (int j = 0; j < i; j++) {
        jump2 = filtered_intersection[queue[j]];
        if (! box_box_intersection(bounding_boxes[jump1], bounding_boxes[jump2]))
          continue;
        for (int k = 0; k < intersect_face[queue[i]].size(); k++) {
          for (int h = 0; h < intersect_face[queue[j]].size(); h++) {
            // todo: move the intersection here
            cut = this->do_intersect(etriangle_eplane,
                                     halfspace[jump1][intersect_face[queue[i]][k]].eplane,
                                     halfspace[jump2][intersect_face[queue[j]][h]].eplane);

            if (!cut) continue;


            cut = is_tpp_on_polyhedra(etriangle_eplane,
                                      halfspace[jump1][intersect_face[queue[i]][k]].eplane,
                                      halfspace[jump2][intersect_face[queue[j]][h]].eplane,
                                      jump1, intersect_face[queue[i]][k]);

            if (!cut) continue;


            cut = is_tpp_on_polyhedra(etriangle_eplane,
                                      halfspace[jump1][intersect_face[queue[i]][k]].eplane,
                                      halfspace[jump2][intersect_face[queue[j]][h]].eplane,
                                      jump2, intersect_face[queue[j]][h]);

            if (!cut) continue;



            inter = Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order(etriangle_eplane,
                                                                                                  halfspace[jump1][intersect_face[queue[i]][k]].eplane,
                                                                                                  halfspace[jump2][intersect_face[queue[j]][h]].eplane,
                                                                                                  idlist, idlistorder, jump1, jump2, check_id);

            if (inter == 1) {


              inter = Implicit_Tri_Facet_Facet_interpoint_Out_Prism_return_local_id_with_face_order_jump_over(etriangle_eplane,
                                                                                                              halfspace[jump1][intersect_face[queue[i]][k]].eplane,

                                                                                                              halfspace[jump2][intersect_face[queue[j]][h]].eplane,
                                                                                                              neighbours, neighbour_facets, neighbour_cover, jump1, jump2, check_id);



              if (inter == 1) {

                return true;
              }
              if (inter == 0) {
                idlistorder.emplace_back(intersect_face[list[check_id]]);
                queue.emplace_back(list[check_id]);
                idlist.emplace_back(filtered_intersection[list[check_id]]);
                coverlist[list[check_id]] = true;
              }
            }
          }
        }
      }
    }

    return false;
  }


  Plane
  get_corner_plane(const Point_3& p0,
                   const Point_3& midp,
                   const Vector_3 &normal,
                   const double distance,
                   const bool robust) const
  {
    Point_3 plane0, plane1, plane2;
    double distance_small = distance * 1;// to be conservative to reduce numerical error, can set the Scalar as 0.999
    Vector_3 direction = CGAL::normalize(p0 - midp);
    plane0 = p0 + direction * distance_small;
    plane1 = plane0 + normal;
    Vector_3 axis =
      //(robust) ? robust_CGAL::cross_product_direction(midp, p0, CGAL::ORIGIN, normal) :
      CGAL::normalize(CGAL::cross_product(direction, normal));
    plane2 = plane0 + axis;

    return Plane(plane0, plane1,plane2);
  }


  // build prisms for a list of triangles. each prism is represented by 7-8 planes, which are represented by 3 points
  void
  halfspace_generation(const std::vector<Point_3> &ver, const std::vector<Vector3i> &faces,
                       std::vector<Prism>& halfspace,
                       std::vector<Iso_cuboid_3>& bounding_boxes, const double &epsilon)
  {
    double tolerance = epsilon / sqrt(3);// the envelope thickness, to be conservative
    double bbox_tolerance = epsilon *(1 + 1e-6);
    Vector_3 AB, AC, BC, normal;
    int de;
    Plane plane;
    std::array<Vector_3, 8> box;


    static const std::array<Vector_3, 8> boxorder = {
      {
        {1, 1, 1},
        {-1, 1, 1},
        {-1, -1, 1},
        {1, -1, 1},
        {1, 1, -1},
        {-1, 1, -1},
        {-1, -1, -1},
        {1, -1, -1},
      } };
    bool use_accurate_cross = false;

    static const int c_face[6][3] = { {0, 1, 2}, {4, 7, 6}, {0, 3, 4}, {1, 0, 4}, {1, 5, 2}, {2, 6, 3} };

    halfspace.resize(faces.size());
    bounding_boxes.resize(faces.size());
    for (int i = 0; i < faces.size(); ++i)
      {
        Bbox_3 bb = ver[faces[i][0]].bbox () + ver[faces[i][1]].bbox() + ver[faces[i][2]].bbox();
        // todo: Add a grow() function to Bbox_3
        bounding_boxes[i] = Iso_cuboid_3(Point_3(bb.xmin()-bbox_tolerance, bb.ymin()-bbox_tolerance, bb.zmin()-bbox_tolerance),
                                         Point_3(bb.xmax()+bbox_tolerance, bb.ymax()+bbox_tolerance, bb.zmax()+bbox_tolerance));

        AB = ver[faces[i][1]] - ver[faces[i][0]];
        AC = ver[faces[i][2]] - ver[faces[i][0]];
        BC = ver[faces[i][2]] - ver[faces[i][1]];

#ifdef TRACE
        de = algorithms::is_triangle_degenerated(ver[faces[i][0]], ver[faces[i][1]], ver[faces[i][2]]);

        if (de == DEGENERATED_POINT)
          {
            //logger().debug("Envelope Triangle Degeneration- Point");
            for (int j = 0; j < 8; j++)
              {
                box[j] = ver[faces[i][0]] + boxorder[j] * tolerance;
              }
            halfspace[i].resize(6);
            for (int j = 0; j < 6; j++) {
              halfspace[i][j][0] = box[c_face[j][0]];
              halfspace[i][j][1] = box[c_face[j][1]];
              halfspace[i][j][2] = box[c_face[j][2]];
            }


            continue;
          }
        if (de == DEGENERATED_SEGMENT)
          {
            //logger().debug("Envelope Triangle Degeneration- Segment");
            Scalar length1 = AB.dot(AB), length2 = AC.dot(AC), length3 = BC.dot(BC);
            if (length1 >= length2 && length1 >= length3)
              {
                algorithms::seg_cube(ver[faces[i][0]], ver[faces[i][1]], tolerance, box);

              }
            if (length2 >= length1 && length2 >= length3)
              {
                algorithms::seg_cube(ver[faces[i][0]], ver[faces[i][2]], tolerance, box);

              }
            if (length3 >= length1 && length3 >= length2)
              {
                algorithms::seg_cube(ver[faces[i][1]], ver[faces[i][2]], tolerance, box);
              }
            halfspace[i].resize(6);
            for (int j = 0; j < 6; j++) {
              halfspace[i][j][0] = box[c_face[j][0]];
              halfspace[i][j][1] = box[c_face[j][1]];
              halfspace[i][j][2] = box[c_face[j][2]];
            }


            continue;
          }
        if (de == NERLY_DEGENERATED)
          {
            //logger().debug("Envelope Triangle Degeneration- Nearly");
            use_accurate_cross = true;

            normal = algorithms::accurate_normal_vector(ver[faces[i][0]], ver[faces[i][1]], ver[faces[i][2]]);

          }
        else
          {
            normal = CGAL::normalize(CGAL::cross_product(AB, AC));
          }
#endif
        normal = CGAL::normalize(CGAL::cross_product(AB, AC)); // remove as soon as #if 1 above

        halfspace[i].reserve(8);
        Vector_3 normaldist = normal * tolerance;
        Vector_3 edgedire, edgenormaldist;
        plane = Plane(ver[faces[i][0]] + normaldist,
                      ver[faces[i][1]] + normaldist,
                      ver[faces[i][2]] + normaldist);
        halfspace[i].emplace_back(plane);// number 0

        plane = Plane(ver[faces[i][0]] - normaldist,
                      ver[faces[i][2]] - normaldist,
                      ver[faces[i][1]] - normaldist);// order: 0, 2, 1
        halfspace[i].emplace_back(plane);// number 1

        int obtuse = CGAL::obtuse_angle(ver[faces[i][0]], ver[faces[i][1]], ver[faces[i][2]]);


        edgedire = CGAL::normalize(AB);
        // if (use_accurate_cross)edgenormaldist = accurate_cross_product_direction(CGAL::ORIGIN, edgedire, CGAL::ORIGIN, normal)*tolerance;
        // else
        edgenormaldist = CGAL::normalize(CGAL::cross_product(edgedire,normal))*tolerance;
        plane = Plane(ver[faces[i][0]] + edgenormaldist,
                      ver[faces[i][1]] + edgenormaldist,
                      ver[faces[i][0]] + edgenormaldist + normal);
        halfspace[i].emplace_back(plane);// number 2


        if (obtuse != 1) {
          plane = get_corner_plane(ver[faces[i][1]], CGAL::midpoint(ver[faces[i][0]], ver[faces[i][2]]) , normal,
                           tolerance, use_accurate_cross);
          halfspace[i].emplace_back(plane);// number 3;

        }

        edgedire = CGAL::normalize(BC);
        // if (use_accurate_cross)edgenormaldist = accurate_cross_product_direction(CGAL::ORIGIN, edgedire, CGAL::ORIGIN, normal)*tolerance;
        // else
        edgenormaldist = CGAL::normalize(CGAL::cross_product(edgedire, normal))*tolerance;

        plane = Plane(ver[faces[i][1]] + edgenormaldist,
                      ver[faces[i][2]] + edgenormaldist,
                      ver[faces[i][1]] + edgenormaldist + normal);
        halfspace[i].emplace_back(plane);// number 4

        if (obtuse != 2) {
          plane = get_corner_plane(ver[faces[i][2]], CGAL::midpoint(ver[faces[i][0]], ver[faces[i][1]]), normal,
                           tolerance,use_accurate_cross);
          halfspace[i].emplace_back(plane);// number 5;

        }

        edgedire = -CGAL::normalize(AC);
        // if (use_accurate_cross)edgenormaldist = accurate_cross_product_direction(CGAL::ORIGIN, edgedire, CGAL::ORIGIN , normal)*tolerance;
        // else
        edgenormaldist = CGAL::normalize(CGAL::cross_product(edgedire, normal))*tolerance;

        plane = Plane(ver[faces[i][2]] + edgenormaldist,
                      ver[faces[i][0]] + edgenormaldist,
                      ver[faces[i][0]] + edgenormaldist + normal);
        halfspace[i].emplace_back(plane);// number 6

        if (obtuse != 0) {
          plane = get_corner_plane(ver[faces[i][0]], CGAL::midpoint(ver[faces[i][1]], ver[faces[i][2]]) , normal,
                           tolerance,use_accurate_cross);
          halfspace[i].emplace_back(plane);// number 7;

        }

#ifdef TRACE
        std::cout << "face "<< i << std::endl;
        for(int j = 0; j < halfspace[i].size(); j++){
          const Plane& p =  halfspace[i][j];
          std::cout << p.ep << " | "  << p.eq << " | "  << p.er << std::endl;
          ePoint_3 pv(ver[faces[i][0]].x(), ver[faces[i][0]].y(),ver[faces[i][0]].z());
          CGAL::Orientation ori = CGAL::orientation(p.ep, p.eq, p.er, pv);
          assert(ori == CGAL::NEGATIVE);
        }
#endif

      }
  }

  void prism_to_off(unsigned int i, std::string fname) const
  {
    std::vector<ePlane_3> eplanes;
    for(int j = 0; j < halfspace[i].size(); j++){
      eplanes.push_back(halfspace[i][j].eplane);
    }
    ePoint_3 origin(env_vertices[env_faces[i][0]].x(), env_vertices[env_faces[i][0]].y(),env_vertices[env_faces[i][0]].z());
    CGAL::Surface_mesh<ePoint_3> esm;
    CGAL::halfspace_intersection_3(eplanes.begin(),eplanes.end(),esm , boost::make_optional(origin));

    CGAL::Surface_mesh<typename CGAL::Exact_predicates_inexact_constructions_kernel::Point_3> sm;
    CGAL::copy_face_graph(esm,sm);
    CGAL::Polygon_mesh_processing::triangulate_faces(sm);
    fname += "_";
    fname += std::to_string(i);
    fname += ".off";
    std::ofstream out(fname.c_str());
    out << sm << std::endl << std::endl;
  }


  // \returns `true` if the query point is inside the envelope
  bool
  operator()(const Point_3& query) const
  {
    std::vector<unsigned int> prismindex;
    tree.all_intersected_primitives(query, std::back_inserter(prismindex));
    if(prismindex.empty()){
      return false;
    }
    ePoint_3 equery(query.x(),query.y(),query.z());
    if(point_out_prism(equery, prismindex, -1)){
      return false;
    }
    return true;
  }


  // \returns `true` if the query segment is inside the envelope
  bool
  operator()(const Point_3& source, const Point_3& target) const
  {
    std::vector<unsigned int> prismindex;
    Segment_3 query(source,target);
    tree.all_intersected_primitives(query, std::back_inserter(prismindex));

    if(segment_out_of_envelope(source, target, prismindex)){
      return false;
    }
    return true;
  }


  // \returns `true` if the query triangle is inside the envelope
  bool
  operator()(const Point_3& t0, const Point_3& t1, const Point_3& t2) const
  {
    std::vector<unsigned int> prismindex;
    Triangle_3 query(t0, t1, t2);
    tree.all_intersected_primitives(query, std::back_inserter(prismindex));
    std::sort(prismindex.begin(), prismindex.end());
#ifdef TRACE
    for(int i=0; i < prismindex.size(); i++){
      std::cout << prismindex[i] << " ";
    }
    std::cout << std::endl;
#endif
    for(int i = 0; i < prismindex.size(); i++){
      prism_to_off(prismindex[i], "prism");
    }
    std::array<Point_3,3> triangle = { t0, t1, t2 };
    if(triangle_out_of_envelope(triangle, prismindex)){
      return false;
    }
    return true;
  }


}; // class Envelope


// `fast` takes an off file and an offset as arguments
//  If called additionally with 3 more vertex indices it performs the envelope test with the triangle
//  Otherwise it tests for all vertex triples forming a non-degenerate trianges
//  and writes the triple in the file inside.txt or outside.txt

int main(int argc, char* argv[])
{
  typedef Kernel::Point_3 Point_3;
  std::vector<Point_3> env_vertices;
  std::vector<Vector3i> env_faces;

  std::ifstream in(argv[1]);

  double eps = std::stod(std::string(argv[2]));

  int ii, ij, ik;

  int query_dim = -1;

  if(argc >= 4){
    query_dim++;
    ii = std::stoi(std::string(argv[3]));
  }

  if(argc >= 5){
    query_dim++;
    ij = std::stoi(std::string(argv[4]));
  }

  if(argc == 6){
    query_dim++;
    ik = std::stoi(std::string(argv[5]));
  }
  std::string off;
  int V, F, E;
  in >> off >> V >> F >> E;
  env_vertices.reserve(V);
  env_faces.reserve(F);

  Kernel::Point_3 p;
  for(int i =0; i < V; i++){
    in >> p;
    env_vertices.push_back(p);
  }
  int three, vi, vj, vk;
  for(int i =0; i < F; ++i){
    in >> three >> vi >> vj >> vk;
    Vector3i f = { vi, vj, vk };
    env_faces.push_back(f);
  }

  CGAL::Timer t;
  t.start();

  //  CGAL::Protect_FPU_rounding<true> pr;

  Envelope<Kernel> envelope(env_vertices, env_faces, eps);


  std::cout << t.time() << " sec." << std::endl;
  t.reset();


  if(query_dim != -1){

    bool bbb;

    if(query_dim == 0){
    Point_3 v0 = env_vertices[ii];
    bbb = envelope(v0);
    } else if (query_dim == 1){
      Point_3 v0 = env_vertices[ii];
      Point_3 v1 = env_vertices[ij];
      bbb = envelope(v0, v1);
    } else {
      Point_3 v0 = env_vertices[ii];
      Point_3 v1 = env_vertices[ij];
      Point_3 v2 = env_vertices[ik];
      {
        std::ofstream query("query.off");
        query << "OFF\n" << "3 1 0\n" << v0 << std::endl << v1 << std::endl << v2 << std::endl << "3 0 1 2" << std::endl;
      }
      bbb = envelope(v0, v1, v2);
    }

    if(bbb){
      std::cout <<  "inside the envelope" << std::endl;
    }else{
      std::cout <<  "outside the envelope" << std::endl;
    }
    return 0;
  }


  std::ofstream inside("inside.txt");
  std::ofstream outside("outside.txt");
  for(int i = 0; i < env_vertices.size(); i+=10){
      for(int j = 0; j < env_vertices.size(); j+= 10){
        for(int k = 0; k < env_vertices.size(); k+=10){
          if( ( i != j) && (i != k) && (j != k)){
            if(! CGAL::collinear(env_vertices[i], env_vertices[j],env_vertices[k])){
              if(envelope(env_vertices[i],  env_vertices[j], env_vertices[k])){
                inside << i << " " << j << " "<< k <<std::endl;
              } else{
                outside << i << " " << j << " "<< k <<std::endl;
              }
            }
          }
        }
      }
  }


  return 0;
}
