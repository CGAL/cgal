#ifndef CGAL_EXACTNESS_H
#define CGAL_EXACTNESS_H
#include <boost/mpl/has_xxx.hpp>
#include <CGAL/tags.h>
namespace CGAL {
namespace internal {
BOOST_MPL_HAS_XXX_TRAIT_DEF(Is_exact)
}
template<class T,bool=internal::has_Is_exact<T>::value> struct Is_exact {
	enum { value=false };
	typedef Tag_false type;
};
template<class T> struct Is_exact<T,true> {
	typedef typename T::Is_exact type;
	enum { value=type::value };
};

}
#endif // CGAL_EXACTNESS_H
