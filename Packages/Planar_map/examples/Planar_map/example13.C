/*! This example demonstrates the technique used to extend the segment
 * kernel-object.
 */
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <iostream>
#include <iterator>
#include <algorithm>

// Kernel is the new kernel, and Kernel_base is the old kernel
template <typename Kernel, typename Kernel_base>
class My_cartesian_base : public Kernel_base::template Base<Kernel>::Type {
public:
  typedef typename Kernel_base::template Base<Kernel>::Type     Old_kernel;
  typedef typename Old_kernel::Segment_2                        Old_segment_2;
  typedef typename Old_kernel::Point_2                          Point_2;
  typedef typename Old_kernel::FT                               FT;

public:
  /*! The new segment type */
  template <class K>
  class New_segment_2 : public K::Old_segment_2 {
  private:
    /*! A field extending the segment */
    int m_data;
    typedef typename K::Point_2 Point_2;
  public:
    /*! Constructors: */
    New_segment_2(int data = 0) : m_data(data) {}
    New_segment_2(const Point_2 & p, const Point_2 & q, int data = 0) :
      Old_segment_2(p, q), m_data(data) {}
    New_segment_2(const Old_segment_2 & seg, int data = 0) :
      Old_segment_2(seg), m_data(data) {}

    /*! \brief obtains the data of the extended segment */
    int get_data() const { return m_data; }
  };
  
  typedef New_segment_2<Kernel>                 Segment_2;

  /*! */
  template <typename Kernel2>
  struct Base { typedef My_cartesian_base<Kernel2, Kernel_base> Type; };
};

/*! The extended Kernel type */
template <typename FT_>
struct Ext_seg_kernel
  : public My_cartesian_base<Ext_seg_kernel<FT_>, CGAL::Cartesian<FT_> >
{};

// The remaining types:
typedef CGAL::Quotient<CGAL::MP_Float>          NT;
typedef Ext_seg_kernel<NT>                      Kernel;
typedef CGAL::Pm_segment_traits_2<Kernel>       Traits;
typedef Traits::Point_2                         Point_2;
typedef Traits::X_monotone_curve_2              X_monotone_curve_2;
typedef CGAL::Pm_default_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>         Planar_map;

/*! Exporter of the extended-curve type */
inline
std::ostream & operator<<(std::ostream & o, const X_monotone_curve_2 & cv)
{
  typedef Kernel::Construct_vertex_2 Construct_vertex_2;
  
  Traits traits;
  Construct_vertex_2 construct_vertex = traits.construct_vertex_2_object();
  const Point_2 & source = construct_vertex(cv, 0);
  const Point_2 & target = construct_vertex(cv, 1);
  o << source << ", " << target << ", " << cv.get_data();
  return o;
}

/*! main */
int main()
{
  // Create an instance of a Planar_map:
  Planar_map pm;
  X_monotone_curve_2 cv[5];

  Point_2 p0(1, 4), p1(5, 7), p2(9, 4), p3(5, 1);

  // Create the curves:
  cv[0] = X_monotone_curve_2(p0, p1, 0);
  cv[1] = X_monotone_curve_2(p1, p2, 1);
  cv[2] = X_monotone_curve_2(p2, p3, 2);
  cv[3] = X_monotone_curve_2(p3, p0, 3);
  cv[4] = X_monotone_curve_2(p0, p2, 4);

  std::cout << "The curves of the map :" << std::endl;
  std::copy(&cv[0], &cv[5],
            std::ostream_iterator<X_monotone_curve_2>(std::cout, "\n"));
  std::cout << std::endl;

  // Insert the curves into the Planar_map:
  std::cout << "Inserting the curves to the map ... ";
  pm.insert(&cv[0], &cv[5]);
  std::cout << ((pm.is_valid()) ? "map valid!" : "map invalid!") << std::endl
            << std::endl;

  // Print the data along the curves:
  Planar_map::Halfedge_const_iterator eit;
  for (eit = pm.halfedges_begin(); eit != pm.halfedges_end(); ++eit, ++eit) {
    const X_monotone_curve_2 & cv = eit->curve();
    std::cout << cv << std::endl;
  }
  
  return 0;
}
