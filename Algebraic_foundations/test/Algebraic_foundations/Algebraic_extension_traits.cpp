#include <CGAL/basic.h>
#include <CGAL/Algebraic_extension_traits.h>
#include <cassert>

int main(){
    typedef CGAL::Algebraic_extension_traits<int> AET;
    
    typedef AET::Type Type;
    BOOST_STATIC_ASSERT((::boost::is_same<int,Type>::value)); 

    typedef AET::Is_extended Is_extended;
    BOOST_STATIC_ASSERT((::boost::is_same<CGAL::Tag_false,Is_extended>::value));
    
    typedef AET::Normalization_factor Normalization_factor; 
    {
        typedef Normalization_factor::argument_type argument_type;
        BOOST_STATIC_ASSERT((::boost::is_same<argument_type,int>::value));
        typedef Normalization_factor::result_type result_type;
        BOOST_STATIC_ASSERT((::boost::is_same<result_type,int>::value));
        Normalization_factor nfac;
        assert(nfac(3)==1);
    }
    typedef AET::Denominator_for_algebraic_integers DFAI; 
    {
        typedef DFAI::argument_type argument_type;
        BOOST_STATIC_ASSERT((::boost::is_same<argument_type,int>::value));
        typedef DFAI::result_type result_type;
        BOOST_STATIC_ASSERT((::boost::is_same<result_type,int>::value));
        DFAI dfai;
        assert(dfai(3)==1);
    }
}
