#include <CGAL/use.h>
#include <CGAL/Real_embeddable_traits.h>
#include <cassert>


#define CGAL_IS_RET_NULL_FUNCTOR(NAME)                                  \
    {                                                                   \
        typedef RET::NAME NAME;                                         \
        CGAL_USE_TYPE(NAME);                                            \
        CGAL_static_assertion(                                            \
                (::boost::is_same<CGAL::Null_functor,NAME>::value));    \
    }

int main(){
    typedef CGAL::Real_embeddable_traits<void> RET;

    typedef RET::Type Type;
    CGAL_USE_TYPE(Type);
    CGAL_static_assertion((::boost::is_same<void,Type>::value));

    typedef RET::Is_real_embeddable Is_real_embeddable;
    CGAL_USE_TYPE(Is_real_embeddable);
    CGAL_static_assertion((::boost::is_same<CGAL::Tag_false,Is_real_embeddable>::value));

    CGAL_IS_RET_NULL_FUNCTOR(Abs);
    CGAL_IS_RET_NULL_FUNCTOR(Sgn);
    CGAL_IS_RET_NULL_FUNCTOR(Is_finite);
    CGAL_IS_RET_NULL_FUNCTOR(Is_positive);
    CGAL_IS_RET_NULL_FUNCTOR(Is_negative);
    CGAL_IS_RET_NULL_FUNCTOR(Is_zero);
    CGAL_IS_RET_NULL_FUNCTOR(Compare);
    CGAL_IS_RET_NULL_FUNCTOR(To_double);
    CGAL_IS_RET_NULL_FUNCTOR(To_interval);
}
