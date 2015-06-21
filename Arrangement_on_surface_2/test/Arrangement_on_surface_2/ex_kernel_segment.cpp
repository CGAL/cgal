/*! This example demonstrates the technique used to extend the segment
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
    int data() const { return m_data; }

    /*! \brief sets the data of the extended point */
    void set_data(int data) { m_data = data; }
  };
  
  template <typename K, typename OldK>
  class New_construct_segment_2
  {
    typedef typename K::RT              RT;
    typedef typename K::Point_2         Point_2;
    typedef typename K::Line_2          Line_2;
    typedef typename K::Segment_2       Segment_2;
    typedef typename Segment_2::Rep     Rep;

  public:
    typedef Segment_2                   result_type;

    // Note : the CGAL::Return_base_tag is really internal CGAL stuff.
    // Unfortunately it is needed for optimizing away copy-constructions,
    // due to current lack of delegating constructors in the C++ standard.

    Rep operator()(CGAL::Return_base_tag, const Point_2 & p, 
                   const Point_2 & q) const
    { return Rep(p, q); }
    Rep operator()(CGAL::Return_base_tag, const Point_2 & p, 
                   const Point_2 & q, int data) const
      { return Rep(p, q, data); }

    // end of hell

    Segment_2 operator()(const Point_2 & p, const Point_2 & q) const
    { return New_segment_2<K>(p, q, 0); }
    Segment_2 operator()(const Point_2 & p, const Point_2 & q, int data) const
    { return New_segment_2<K>(p, q, data); }
  };

  typedef New_segment_2<Kernel>                         Segment_2;
  typedef New_construct_segment_2<Kernel, Old_kernel>   Construct_segment_2;

  Construct_segment_2 construct_segment_2_object() const
  { return Construct_segment_2(); }
  
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

typedef CGAL::Exact_rational                            Number_type;
typedef Ext_seg_kernel<Number_type>                     Kernel;
typedef CGAL::Arr_non_caching_segment_traits_2<Kernel>  Traits;
typedef Traits::Point_2                                 Point_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits>                     Arr;

/*! Exporter of the extended-curve type */
inline
std::ostream & operator<<(std::ostream & o, const X_monotone_curve_2 & cv)
{
  typedef Kernel::Construct_vertex_2 Construct_vertex_2;
  
  Traits traits;
  Construct_vertex_2 construct_vertex = traits.construct_vertex_2_object();
  const Point_2 & source = construct_vertex(cv, 0);
  const Point_2 & target = construct_vertex(cv, 1);
  o << source << ", " << target << ", " << cv.data();
  return o;
}

/*! main */
int main()
{
  // Create an instance of an arrangement:
  Arr arr;
  X_monotone_curve_2 cv[5];

  Point_2 p0(1, 4), p1(5, 7), p2(9, 4), p3(5, 1);

  // Create the curves:
  cv[0] = X_monotone_curve_2(p0, p1);
  cv[1] = X_monotone_curve_2(p1, p2);
  cv[2] = X_monotone_curve_2(p2, p3);
  cv[3] = X_monotone_curve_2(p3, p0);
  cv[4] = X_monotone_curve_2(p0, p2);
  cv[0].set_data(0);
  cv[1].set_data(1);
  cv[2].set_data(2);
  cv[3].set_data(3);
  cv[4].set_data(4);

  std::cout << "The curves of the map :" << std::endl;
  std::copy(&cv[0], &cv[5],
            std::ostream_iterator<X_monotone_curve_2>(std::cout, "\n"));
  std::cout << std::endl;

  // Insert the curves into the arrangement:
  std::cout << "Inserting the curves to the map ... ";
  insert(arr, &cv[0], &cv[5]);
  std::cout << ((arr.is_valid()) ? "map valid!" : "map invalid!") << std::endl
            << std::endl;

  // Print the data along the curves:
  Arr::Halfedge_const_iterator eit;
  for (eit = arr.halfedges_begin(); eit != arr.halfedges_end(); ++eit, ++eit) {
    const X_monotone_curve_2 & cv = eit->curve();
    std::cout << cv << std::endl;
  }
  
  return 0;
}
