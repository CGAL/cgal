#include <CGAL/use.h>
#include <CGAL/Algebraic_extension_traits.h>
#include <CGAL/Sqrt_extension.h>
#include <cassert>

int main(){
    {
        typedef CGAL::Algebraic_extension_traits<int> AET;

        typedef AET::Type Type;
        CGAL_USE_TYPE(Type);
        CGAL_static_assertion((::boost::is_same<int,Type>::value));

        typedef AET::Is_extended Is_extended;
        CGAL_USE_TYPE(Is_extended);
        CGAL_static_assertion(
                (::boost::is_same<CGAL::Tag_false,Is_extended>::value));

        typedef AET::Normalization_factor Normalization_factor;
        {
            typedef Normalization_factor::argument_type argument_type;
            CGAL_USE_TYPE(argument_type);
            CGAL_static_assertion((::boost::is_same<argument_type,int>::value));
            typedef Normalization_factor::result_type result_type;
            CGAL_USE_TYPE(result_type);
            CGAL_static_assertion((::boost::is_same<result_type,int>::value));
            Normalization_factor nfac;
            assert(nfac(3)==1);
        }
        typedef AET::Denominator_for_algebraic_integers DFAI;
        {
            typedef DFAI::argument_type argument_type;
            CGAL_USE_TYPE(argument_type);
            CGAL_static_assertion((::boost::is_same<argument_type,int>::value));
            typedef DFAI::result_type result_type;
            CGAL_USE_TYPE(result_type);
            CGAL_static_assertion((::boost::is_same<result_type,int>::value));
            DFAI dfai;
            assert(dfai(3)==1);
        }
    }
    {
        typedef CGAL::Sqrt_extension<int,int> EXT;
        typedef CGAL::Algebraic_extension_traits<EXT> AET;

        typedef AET::Type Type;
        CGAL_USE_TYPE(Type);
        CGAL_static_assertion((::boost::is_same<EXT,Type>::value));

        typedef AET::Is_extended Is_extended;
        CGAL_USE_TYPE(Is_extended);
        CGAL_static_assertion(
                (::boost::is_same<CGAL::Tag_true,Is_extended>::value));

        typedef AET::Normalization_factor Normalization_factor;
        {
            typedef Normalization_factor::argument_type argument_type;
            CGAL_USE_TYPE(argument_type);
            CGAL_static_assertion((::boost::is_same<argument_type,EXT>::value));
            typedef Normalization_factor::result_type result_type;
            CGAL_USE_TYPE(result_type);
            CGAL_static_assertion((::boost::is_same<result_type,EXT>::value));
            Normalization_factor nfac;
            assert(nfac(EXT(3))==1);
            assert(nfac(EXT(3,0,5))==1);
            assert(nfac(EXT(3,1,5))==EXT(3,-1,5));
        }
        typedef AET::Denominator_for_algebraic_integers DFAI;
        {
            typedef DFAI::argument_type argument_type;
            CGAL_USE_TYPE(argument_type);
            CGAL_static_assertion((::boost::is_same<argument_type,EXT>::value));
            typedef DFAI::result_type result_type;
            CGAL_USE_TYPE(result_type);
            CGAL_static_assertion((::boost::is_same<result_type,EXT>::value));
            DFAI dfai;
            assert(dfai(EXT(3))==1);
            assert(dfai(EXT(3,0,5))==1);
            assert(dfai(EXT(3,1,5))==20);
        }
    }
}
