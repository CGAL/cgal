#ifndef CGAL_KERNELD_TYPES_ISO_BOX_H
#define CGAL_KERNELD_TYPES_ISO_BOX_H
#include <utility>
#include <CGAL/basic.h>
#include <CGAL/functor_tags.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/transforming_pair_iterator.h>
namespace CGAL {
template <class R_> class Iso_box {
	typedef typename Get_type<R_, FT_tag>::type FT_;
	typedef typename Get_type<R_, Point_tag>::type	Point_;
	typedef std::pair<Point_,Point_> Data_;
	Data_ data;
	public:
	Iso_box(){}
	Iso_box(Point_ const&a, Point_ const&b): data(a,b) {}
	Point_ min()const{
	  return data.first;
	}
	Point_ max()const{
	  return data.second;
	}
};
namespace CartesianDKernelFunctors {
  template <class R_> struct Construct_iso_box : Store_kernel<R_> {
    CGAL_FUNCTOR_INIT_STORE(Construct_iso_box)
    typedef typename Get_type<R_, Iso_box_tag>::type	result_type;
    typedef typename Get_type<R_, RT_tag>::type RT;
    typedef typename Get_type<R_, Point_tag>::type	Point;
    typedef typename Get_type<R_, Vector_tag>::type	Vector;
    typedef typename Get_functor<R_, Construct_ttag<Point_tag> >::type Cp_;
    typedef typename Get_functor<R_, Construct_ttag<Point_cartesian_const_iterator_tag> >::type Ci_;
    result_type operator()(Point const&a, Point const&b)const{
      Cp_ cp(this->kernel());
      Ci_ ci(this->kernel());
      return result_type(cp(
	  make_transforming_pair_iterator(ci(a,Begin_tag()), ci(b,Begin_tag()), Min<RT>()),
	  make_transforming_pair_iterator(ci(a,End_tag()), ci(b,End_tag()), Min<RT>())),
      cp(
	  make_transforming_pair_iterator(ci(a,Begin_tag()), ci(b,Begin_tag()), Max<RT>()),
	  make_transforming_pair_iterator(ci(a,End_tag()), ci(b,End_tag()), Max<RT>())));
    }
  };
}

} // namespace CGAL

#endif // CGAL_KERNELD_TYPES_ISO_BOX_H
