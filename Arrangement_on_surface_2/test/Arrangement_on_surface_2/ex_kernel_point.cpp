/*! This example demonstrates the technique used to extend the point
 * kernel-object.
 */
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Arr_non_caching_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <iostream>
#include <iterator>
#include <algorithm>

// Kernel is the new kernel, and Kernel_base is the old kernel
template <typename Kernel, typename Kernel_base>
class My_cartesian_base : public Kernel_base::template Base<Kernel>::Type {
public:
  typedef typename Kernel_base::template Base<Kernel>::Type     Old_kernel;
  typedef typename Old_kernel::Point_2                          Old_point_2;
  typedef typename Old_kernel::FT                               FT;
  typedef typename Old_kernel::RT                               RT;
  typedef typename Old_kernel::Segment_2                        Segment_2;

public:
  /*! The new point type */
  class New_point_2 : public Old_point_2 {
  private:
    /*! A field extending the point */
    int m_data;

  public:
    /*! Constructors */
    New_point_2(int data = 0) : m_data(data) {}
    New_point_2(const CGAL::Origin & origin, int data = 0) :
      Old_point_2(origin), m_data(data) {}

    New_point_2(const RT & hx, const RT & hy, int data = 0) :
      Old_point_2(hx, hy), m_data(data) {}

    /*! \brief obtains the data of the extended point */
    int data() const { return m_data; }

    /*! \brief sets the data of the extended point */
    void set_data(int data) { m_data = data; }
  };

  template <typename K, typename OldK> class New_construct_point_2 {
    typedef typename K::RT              RT;
    typedef typename K::Point_2         Point_2;
    typedef typename K::Line_2          Line_2;
    typedef typename Point_2::Rep       Rep;


  public:
    typedef Point_2                     result_type;

    // Note : the CGAL::Return_base_tag is really internal CGAL stuff.
    // Unfortunately it is needed for optimizing away copy-constructions,
    // due to current lack of delegating constructors in the C++ standard.
    Rep operator()(CGAL::Return_base_tag, CGAL::Origin o) const
    { return Rep(o); }

    Rep operator()(CGAL::Return_base_tag, const RT& x, const RT& y) const
    { return Rep(x, y); }

    Rep operator()(CGAL::Return_base_tag, const RT& x, const RT& y,
                   const RT& w) const
    { return Rep(x, y, w); }

    // End of hell

    Point_2 operator()(CGAL::Origin o) const { return New_point_2(0, 0, 0); }

    Point_2 operator()(const RT & x, const RT & y) const
    { return New_point_2(x, y, 0); }

    Point_2 operator()(const Line_2 & l) const
    {
      typename OldK::Construct_point_2 base_operator;
      Point_2 p = base_operator(l);
      return p;
    }

    Point_2 operator()(const Line_2 & l, int i) const
    {
      typename OldK::Construct_point_2 base_operator;
      return base_operator(l, i);
    }

    // We need this one, as such a functor is in the Filtered_kernel
    Point_2 operator()(const RT & x, const RT & y, const RT & w) const
    {
      if(w != 1) {
        return New_point_2(x/w, y/w, 0);
      } else {
        return New_point_2(x,y, 0);
      }
    }
  };

  typedef New_point_2                                   Point_2;
  typedef New_construct_point_2<Kernel, Old_kernel>     Construct_point_2;

  Construct_point_2 construct_point_2_object() const
  { return Construct_point_2(); }

  /*! */
  template <typename Kernel2>
  struct Base { typedef My_cartesian_base<Kernel2, Kernel_base> Type; };
};

/*! The extended Kernel type */
template <typename FT_>
struct Ext_seg_kernel :
  public CGAL::Type_equality_wrapper<My_cartesian_base<Ext_seg_kernel<FT_>,
                                                       CGAL::Cartesian<FT_> >,
                                     Ext_seg_kernel<FT_> >
{};

// The remaining types:
// leda_rational, or Gmpq, or Quotient<MP_float>
typedef CGAL::Exact_rational                            Number_type;
typedef Ext_seg_kernel<Number_type>                     Kernel;
typedef CGAL::Arr_non_caching_segment_traits_2<Kernel>  Traits;
typedef Traits::Point_2                                 Point_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits>                     Arr;

/*! Exporter of the extended-point type */
inline
std::ostream & operator<<(std::ostream & o, const Point_2 & p)
{
  o << p.x() << ", " << p.y() << ", " << p.data();
  return o;
}

/*! main */
int main()
{
  // Create an instance of a Planar_map:
  Arr arr;
  X_monotone_curve_2 cv[5];

  Point_2 p0(1, 4), p1(5, 7), p2(9, 4), p3(5, 1);
  p0.set_data(0);
  p1.set_data(1);
  p2.set_data(2);
  p3.set_data(3);

  // Create the curves:
  cv[0] = X_monotone_curve_2(p0, p1);
  cv[1] = X_monotone_curve_2(p1, p2);
  cv[2] = X_monotone_curve_2(p2, p3);
  cv[3] = X_monotone_curve_2(p3, p0);
  cv[4] = X_monotone_curve_2(p0, p2);

  std::cout << "The curves of the map :" << std::endl;
  std::copy(&cv[0], &cv[5],
            std::ostream_iterator<X_monotone_curve_2>(std::cout, "\n"));
  std::cout << std::endl;

  // Insert the curves into the Planar_map:
  std::cout << "Inserting the curves to the map ... ";
  insert(arr, &cv[0], &cv[5]);
  std::cout << ((arr.is_valid()) ? "map valid!" : "map invalid!") << std::endl
            << std::endl;

  // Print the curves:
  Arr::Halfedge_const_iterator eit;
  for (eit = arr.halfedges_begin(); eit != arr.halfedges_end(); ++eit, ++eit) {
    const X_monotone_curve_2 & cv = eit->curve();
    std::cout << cv << std::endl;
  }

  return 0;
}
