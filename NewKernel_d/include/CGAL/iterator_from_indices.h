#ifndef CGAL_ITERATOR_FROM_INDICES_H
#define CGAL_ITERATOR_FROM_INDICES_H
#include <boost/iterator/iterator_facade.hpp>
namespace CGAL {
template <class Ref_>
struct Default_coordinate_access {
	typedef Ref_ result_type;
	template<class T> Ref_ operator()(T const& t, int i)const{
		return t[i];
	}
};

//TODO: default type for Value_: typename same_cv<Container_,typename remove_cv<Container_>::type::value_type>::type
template <class Container_, class Value_, class Ref_=
#ifdef CGAL_CXX0X
	decltype(std::declval<Container_>()[0])
#else
	Value_&
#endif
	, class Coord_access = Default_coordinate_access<Ref_>
	>
class Iterator_from_indices
: public boost::iterator_facade<Iterator_from_indices<Container_,Value_,Ref_,Coord_access>,
	Value_, std::bidirectional_iterator_tag, Ref_>
{
	friend class boost::iterator_core_access;
	//FIXME: use int to save space
	//FIXME: use a signed type
	typedef std::size_t index_t;
	Container_& cont;
	index_t index;
	void increment(){ ++index; }
	void decrement(){ --index; }
	void advance(std::ptrdiff_t n){ index+=n; }
	ptrdiff_t distance_to(Iterator_from_indices const& other)const{
		return static_cast<ptrdiff_t>(other.index)
			-static_cast<ptrdiff_t>(index);
	}
	bool equal(Iterator_from_indices const& other)const{
		return index==other.index;
	}
	Ref_ dereference()const{
		//FIXME: use the functor properly
		Coord_access ca;
		return ca(cont,index);
	}
	public:
	Iterator_from_indices(Container_& cont_,std::size_t n)
		: cont(cont_), index(n) {}
};
}
#endif // CGAL_ITERATOR_FROM_INDICES_H
