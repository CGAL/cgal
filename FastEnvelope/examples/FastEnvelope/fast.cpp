
// #define CGAL_PROFILE

#include <Eigen/Dense>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_primitive.h>
#include <boost/iterator/counting_iterator.hpp>
#include <fstream>
#include <string>

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

typedef CGAL::Simple_cartesian<double> Kernel;



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


  std::vector<Point_3> env_vertices;
  std::vector<Vector3i> env_faces;


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
    cid.clear();
    std::vector<bool> cut;
    const Prism& prism = halfspace[cindex];

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


    // todo: Use a lazy plane
    for (int i = 0; i < cutp.size(); i++){
      const Plane& plane_i = prism[cutp[i]];

      CGAL::cpp11::result_of<eIntersect_3(eLine_3, ePlane_3)>::type
        result = CGAL::intersection(line, plane_i.eplane);
      if(! result){
        std::cout <<  "there must be an intersection" << std::endl;
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
    if(! result){
      std::cout <<  "there must be an intersection" << std::endl;
    }

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



  bool segment_out_of_envelope(const Point_3& source, const Point_3& target,
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




  Plane
  get_corner_plane(const Point_3& p0, const Point_3& midp, const Vector_3 &normal, const double distance,
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
  void halfspace_generation(const std::vector<Point_3> &ver, const std::vector<Vector3i> &faces,
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

#if 0
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
        /*
        for(int j = 0; j < halfspace[i].size(); j++){
          const Plane& p =  halfspace[i][j];
          CGAL::Orientation ori = orientation(p.ep, p.eq, p.er, ver[faces[i][0]]);
          assert(ori == CGAL::NEGATIVE);
        }
        */


      }
  }

  // \returns `true` if the query point is inside the envelope
  bool operator()(const Point_3& query) const
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
  bool operator()(const Point_3& source, const Point_3& target) const
  {
    std::vector<unsigned int> prismindex;
    Segment_3 query(source,target);
    tree.all_intersected_primitives(query, std::back_inserter(prismindex));

    if(segment_out_of_envelope(source, target, prismindex)){
      return false;
    }
    return true;
  }

#ifdef TRIANGLE_OUT_OF_ENVELOPE
  // \returns `true` if the query segment is inside the envelope
  bool operator()(const Point_3& t0, const Point_3& t1, const Point_3& t2) const
  {
    std::vector<unsigned int> prismindex;
    Triangle_3 query(t0, t1, t2);
    tree.all_intersected_primitives(query, std::back_inserter(prismindex));
    if(triangle_out_of_envelope(source,target, prismindex)){
      return false;
    }
    return true;
  }
#endif

}; // class Envelope




int main(int argc, char* argv[])
{
  typedef Kernel::Point_3 Point_3;
  std::vector<Point_3> env_vertices;
  std::vector<Vector3i> env_faces;


  std::ifstream in(argv[1]);

  double eps = std::stod(std::string(argv[2]));

  std::string off;
  int V, F, E;
  in >> off >> V >> F >> E;
  env_vertices.reserve(V);
  env_faces.reserve(F);

  std::ofstream vert("vertices.txt");
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

  Point_3 v0 = env_vertices[0];
  Point_3 v44 = env_vertices[44];
  bool b0_44 = envelope(v0,v44);

  std::ofstream out("inside.polylines.txt");
  int count = 0;
  for(int i = 0; i < env_vertices.size(); i++){
    Point_3 vi = env_vertices[i];
    if(envelope(v0,vi)){
      count++;
      out << "2 " << v0 << " " << vi << std::endl;
    }
  }
  std::cout << count << " inside in " << t.time() << " sec." << std::endl;
  return 0;
}
