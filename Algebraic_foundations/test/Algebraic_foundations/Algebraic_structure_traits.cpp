#include <CGAL/use.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <cassert>


#define CGAL_IS_AST_NULL_FUNCTOR(NAME)                                  \
    {                                                                   \
        typedef AST::NAME NAME;                                         \
        CGAL_USE_TYPE(NAME);                                             \
        static_assert(                                            \
                (::std::is_same<CGAL::Null_functor,NAME>::value));    \
    }

int main(){
    typedef CGAL::Algebraic_structure_traits<void> AST;

    typedef AST::Type Type;
    CGAL_USE_TYPE(Type);
    static_assert(::std::is_same<void,Type>::value);

    typedef AST::Algebraic_category Algebraic_category;
    CGAL_USE_TYPE(Algebraic_category);
    static_assert(
        (::std::is_same<CGAL::Null_tag,Algebraic_category>::value));

    typedef AST::Is_exact Is_exact;
    CGAL_USE_TYPE(Is_exact);
    static_assert(::std::is_same<CGAL::Null_tag,Is_exact>::value);
    typedef AST::Is_numerical_sensitive Is_sensitive;
    CGAL_USE_TYPE(Is_sensitive);
    static_assert(::std::is_same<CGAL::Null_tag,Is_sensitive>::value);

    CGAL_IS_AST_NULL_FUNCTOR (  Simplify);
    CGAL_IS_AST_NULL_FUNCTOR (  Unit_part);
    CGAL_IS_AST_NULL_FUNCTOR (  Integral_division);
    CGAL_IS_AST_NULL_FUNCTOR (  Divides);
    CGAL_IS_AST_NULL_FUNCTOR (  Is_square);
    CGAL_IS_AST_NULL_FUNCTOR (  Gcd);
    CGAL_IS_AST_NULL_FUNCTOR (  Div_mod);
    CGAL_IS_AST_NULL_FUNCTOR (  Div);
    CGAL_IS_AST_NULL_FUNCTOR (  Mod);
    CGAL_IS_AST_NULL_FUNCTOR (  Square);
    CGAL_IS_AST_NULL_FUNCTOR (  Is_zero);
    CGAL_IS_AST_NULL_FUNCTOR (  Is_one);
    CGAL_IS_AST_NULL_FUNCTOR (  Sqrt);
    CGAL_IS_AST_NULL_FUNCTOR (  Kth_root);
    CGAL_IS_AST_NULL_FUNCTOR (  Root_of);
}
