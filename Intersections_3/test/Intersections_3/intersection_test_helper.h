#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Tetrahedron_3.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Random.h>

#include "create_bbox_mesh.h"

#include <cassert>
#include <iostream>
#include <vector>

const double epsilon = 0.001;

struct randomint
{
  randomint() ;
  int get() const { return sequence[cur]; }
  int next() { cur = (cur+1)%11; return get();}

private:
  int sequence[11];
  int cur;
};

inline randomint::randomint()
{
  cur = 0;
  sequence[0] = 19;
  sequence[1] = 5;
  sequence[2] = 17;
  sequence[3] = 13;
  sequence[4] = 29;
  sequence[5] = 2;
  sequence[6] = 23;
  sequence[7] = 31;
  sequence[8] = 3;
  sequence[9] = 37;
  sequence[10] = 11;
}

randomint ri;

inline double to_nt(int d)
{
  return double(d);
}

template <typename K>
struct Intersection_3_tester
{
  typedef CGAL::Bbox_3              Bbox;

  typedef CGAL::Iso_cuboid_3<K>     Cub;
  typedef CGAL::Line_3<K>           L;
  typedef CGAL::Point_3<K>          P;
  typedef CGAL::Plane_3<K>          Pl;
  typedef CGAL::Ray_3<K>            R;
  typedef CGAL::Segment_3<K>        S;
  typedef CGAL::Sphere_3<K>         Sph;
  typedef CGAL::Tetrahedron_3<K>    Tet;
  typedef CGAL::Triangle_3<K>       Tr;

  typedef std::vector<P>            Pol;

protected:
  CGAL::Random& r;

  bool has_exact_p;
  bool has_exact_c;

  int m = 0, M = 100;

  bool verbose = false;

public:
  Intersection_3_tester(CGAL::Random& r,
                        const bool has_exact_p = false,
                        const bool has_exact_c = false)
    : r(r), has_exact_p(has_exact_p), has_exact_c(has_exact_c)
  { }

public:
  Pl pl(int a, int b, int c, int d)
  {
    int w = ri.next();
    return Pl(to_nt(a*w), to_nt(b*w), to_nt(c*w), to_nt(d*w));
  }

  P p(int x, int y, int z)
  {
    int w = ri.next();
    return P(to_nt(x*w), to_nt(y*w), to_nt(z*w), to_nt(w));
  }

  P random_point()
  {
    return p(r.get_int(m, M), r.get_int(m, M), r.get_int(m, M));
  }

public:
  template <typename Type>
  bool approx_equal_nt(const Type &t1, const Type &t2)
  {
    if(t1 == t2)
      return true;

    if(has_exact_c)
    {
      std::cerr << " Comparison failed between : " << t1 << " and " << t2 << "\n";
      return false;
    }

    if(CGAL::abs(t1 - t2) / (CGAL::max)(CGAL::abs(t1), CGAL::abs(t2)) < epsilon)
      return true;

    std::cerr << " Approximate comparison failed between : " << t1 << " and " << t2 << "\n";

    return false;
  }

  // we need approx equal to check approx kernels, but maybe we should only test with exact kernels
  // (approx kernels were useful before, when the text output was checked by diff ?)
  // idea : test containment with intervals ?  or use some "epsilon double"?
  // I need to convert the text output to exact rationals in the source...
  // Well, for now the current scheme works.
  template <typename Type>
  bool approx_equal(const Type& t1, const Type& t2)
  {
    if(t1 != t2)
    {
      if(verbose)
        std::cerr << "Failed comparison of " << typeid(Type).name() << ": " << t1 << " is NOT " << t2 << std::endl;

      return false;
    }

    return true;
  }

//  bool approx_equal(const P& p, const P& q)
//  {
//    return approx_equal_nt(p.x(), q.x()) &&
//           approx_equal_nt(p.y(), q.y()) &&
//           approx_equal_nt(p.z(), q.z());
//  }

  bool approx_equal(const S& p, const S& q)
  {
    // allow opposite orientation
    return (approx_equal(p.source(), q.source()) && approx_equal(p.target(), q.target())) ||
           (approx_equal(p.source(), q.target()) && approx_equal(p.target(), q.source()));
  }

  bool approx_equal(const Pol& p, const Pol& q)
  {
    if(p.size() != q.size())
      return false;

    for(typename Pol::const_iterator itp = p.begin(), itq = q.begin(); itp != p.end(); ++itp, ++itq)
      if(!approx_equal(*itp, *itq))
        return false;

    return true;
  }

public:
  template <typename O1, typename O2>
  void check_do_intersect(const O1& o1, const O2& o2)
  {
    if(verbose)
    {
      std::cout << "\nExpecting intersection between " << typeid(O1).name() << " = " << o1
                << "\nand " << typeid(O2).name() << " = " << o2 << std::endl;
    }

    const bool res12 = CGAL::do_intersect(o1, o2);
    const bool res21 = CGAL::do_intersect(o2, o1);

    if(has_exact_p)
    {
      assert(res12);
      assert(res21);
    }
    else
    {
      CGAL_warning(res12);
      CGAL_warning(res21);
    }
  }

  template <typename Res, typename O1, typename O2>
  void check_intersection(const O1& o1, const O2& o2)
  {
    check_do_intersect(o1, o2);

    const auto res12 = CGAL::intersection(o1, o2);
    const auto res21 = CGAL::intersection(o2, o1);

    Res tmp;
    if(has_exact_p)
    {
      assert(CGAL::assign(tmp, res12));
      assert(CGAL::assign(tmp, res21));
    }
    else
    {
      CGAL_warning(CGAL::assign(tmp, res12));
      CGAL_warning(CGAL::assign(tmp, res21));
    }
  }

  template <typename Res, typename O1, typename O2>
  void check_intersection(const O1& o1, const O2& o2,
                          const Res& result,
                          bool do_opposite = true)
  {
    if(verbose)
    {
      std::cout << "\nCheck intersection between " << typeid(O1).name() << " = " << o1
                << "\nand " << typeid(O2).name() << " = " << o2 << std::endl;
      std::cout << "Expecting result " << typeid(Res).name() << " = " << result << std::endl;
    }

    const bool res12 = CGAL::do_intersect(o1, o2);
    if(has_exact_p)
      assert(res12);
    else
      CGAL_warning(res12);

    const auto ires12 = CGAL::intersection(o1, o2);

    Res tmp;
    if(has_exact_c)
    {
      assert(CGAL::assign(tmp, ires12));
      assert(approx_equal(tmp, result));
    }
    else
    {
      if(CGAL::assign(tmp, ires12))
        CGAL_warning(approx_equal(tmp, result));
      else
        CGAL_warning_msg(false, "Expected an intersection, but it was not found!");
    }

    if(do_opposite)
      check_intersection(o2, o1, result, false);
  }

  template <typename Res, typename O1, typename O2, typename O3, typename O4>
  void check_intersection(const O1& o1, const O2& o2,
                          const O3& o3, const O4& o4)
  {
    if(verbose)
    {
      std::cout << "\nCheck intersection between " << typeid(O1).name() << " = " << o1
                << "\nand " << typeid(O2).name() << " = " << o2 << std::endl;
      std::cout << "Expecting result to be the same as intersection between " << typeid(O3).name() << " = " << o3
                << "\nand " << typeid(O4).name() << " = " << o4 << std::endl;
    }

    const bool res12 = CGAL::do_intersect(o1, o2);
    const bool res34 = CGAL::do_intersect(o3, o4);

    if(has_exact_p)
      assert(res12 == res34);
    else
      CGAL_warning(res12 == res34);

    if(res12 && res34)
    {
      const auto ires12 = CGAL::intersection(o1, o2);
      const auto ires34 = CGAL::intersection(o3, o4);

      Res tmp12, tmp34;
      if(has_exact_c)
      {
        assert(CGAL::assign(tmp12, ires12));
        assert(CGAL::assign(tmp34, ires34));
        assert(tmp12 == tmp34);
      }
      else
      {
        if(CGAL::assign(tmp12, ires12) && CGAL::assign(tmp34, ires34))
        {
          CGAL_warning(tmp12 == tmp34);
        }
      }
    }
  }

  template <typename OV>
  struct Variant_visitor
  {
    Variant_visitor(const OV& other_variant)
      : ov(other_variant), equal(false)
    {
      assert(ov);
    }

    template<typename T>
    bool compare_to_other_variant(const T& t) const
    {
      if(auto* r = std::get_if<T>(&*ov))
      {
        return (t == *r);
      }

      return false;
    }

    template <typename T>
    void operator()(const T& value) const
    {
      if(equal)
        return;

      if(compare_to_other_variant(value))
        equal = true;
    }

    const OV& ov;
    mutable bool equal;
  };

  template <typename O1, typename O2, typename O3, typename O4>
  void check_intersection(const O1& o1, const O2& o2,
                          const O3& o3, const O4& o4)
  {
    if(verbose)
    {
      std::cout << "\nCheck intersection between " << typeid(O1).name() << " = " << o1
                << "\nand " << typeid(O2).name() << " = " << o2 << std::endl;
      std::cout << "Expecting result to be the same as intersection between " << typeid(O3).name() << " = " << o3
                << "\nand " << typeid(O4).name() << " = " << o4 << std::endl;
    }

    const bool res12 = CGAL::do_intersect(o1, o2);
    const bool res34 = CGAL::do_intersect(o3, o4);

    if(has_exact_p)
      assert(res12 == res34);
    else
      CGAL_warning(res12 == res34);

    if(res12 && res34)
    {
      const auto ires12 = CGAL::intersection(o1, o2);
      const auto ires34 = CGAL::intersection(o3, o4);

      if(has_exact_c)
      {
        assert(ires12 && ires34);

        Variant_visitor<decltype(ires12)> vis(ires12);
        std::visit(vis, *ires34);
        assert(vis.equal);
      }
    }
  }

public:
  template <typename O1, typename O2>
  void check_do_not_intersect(const O1& o1, const O2& o2)
  {
    const bool res12 = CGAL::do_intersect(o1, o2);
    const bool res21 = CGAL::do_intersect(o2, o1);

    if(has_exact_p)
    {
      assert(!res12);
      assert(!res21);
    }
    else
    {
      CGAL_warning(!res12);
      CGAL_warning(!res21);
    }
  }

  template <typename O1, typename O2, typename O3>
  void check_do_not_intersect(const O1& o1, const O2& o2, const O3& o3)
  {
    auto res123 = CGAL::do_intersect(o1, o2, o3);
    auto res132 = CGAL::do_intersect(o1, o3, o2);
    auto res213 = CGAL::do_intersect(o2, o1, o3);
    auto res231 = CGAL::do_intersect(o2, o3, o1);
    auto res312 = CGAL::do_intersect(o3, o1, o2);
    auto res321 = CGAL::do_intersect(o3, o2, o1);
    if(has_exact_p)
      assert(!res123 && !res132 && !res213 && !res231 && !res312 && !res321);
    else
      CGAL_warning(!res123 && !res132 && !res213 && !res231 && !res312 && !res321);
  }

  template <typename O1, typename O2>
  void check_no_intersection(const O1& o1, const O2& o2)
  {
    check_do_not_intersect(o1, o2);

    auto res12 = CGAL::intersection(o1, o2);
    auto res21 = CGAL::intersection(o2, o1);
    if(has_exact_p)
    {
      assert(!res12);
      assert(!res21);
    }
    else
    {
      CGAL_warning(!res12);
      CGAL_warning(!res21);
    }
  }

  template <typename O1, typename O2, typename O3>
  void check_no_intersection(const O1& o1, const O2& o2, const O3& o3)
  {
    check_do_not_intersect(o1, o2, o3);

    auto res123 = CGAL::intersection(o1, o2, o3);
    auto res132 = CGAL::intersection(o1, o3, o2);
    auto res213 = CGAL::intersection(o2, o1, o3);
    auto res231 = CGAL::intersection(o2, o3, o1);
    auto res312 = CGAL::intersection(o3, o1, o2);
    auto res321 = CGAL::intersection(o3, o2, o1);
    if(has_exact_p)
      assert(!res123 && !res132 && !res213 && !res231 && !res312 && !res321);
    else
      CGAL_warning(!res123 && !res132 && !res213 && !res231 && !res312 && !res321);
  }
};
