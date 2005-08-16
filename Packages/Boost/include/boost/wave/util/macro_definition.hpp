/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library

    http://www.boost.org/

    Copyright (c) 2001-2005 Hartmut Kaiser. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(MACRO_DEFINITION_HPP_D68A639E_2DA5_4E9C_8ACD_CFE6B903831E_INCLUDED)
#define MACRO_DEFINITION_HPP_D68A639E_2DA5_4E9C_8ACD_CFE6B903831E_INCLUDED

#include <vector>
#include <list>

#include <boost/wave/wave_config.hpp>
#include <boost/wave/token_ids.hpp>

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave {
namespace util {

///////////////////////////////////////////////////////////////////////////////
//
//  macro_definition
//
//      This class containes all infos for a defined macro. 
//
///////////////////////////////////////////////////////////////////////////////
template <typename TokenT, typename ContainerT>
struct macro_definition {

    typedef std::vector<TokenT> parameter_container_t;
    typedef ContainerT          definition_container_t;
    
    typedef typename parameter_container_t::const_iterator 
        const_parameter_iterator_t;
    typedef typename definition_container_t::const_iterator 
        const_definition_iterator_t;

    macro_definition(TokenT const &token_, bool has_parameters, 
            bool is_predefined_, long uid_)
    :   macroname(token_), uid(uid_), is_functionlike(has_parameters), 
        replaced_parameters(false), is_available_for_replacement(true),
        is_predefined(is_predefined_)
#if BOOST_WAVE_SUPPORT_VARIADICS_PLACEMARKERS != 0
        , has_ellipsis(false)
#endif
    {
    }
    // generated copy constructor
    // generated destructor
    // generated assignment operator

    // Replace all occurrences of the parameters throughout the macrodefinition
    // with special parameter tokens to simplify later macro replacement.
    // Additionally mark all occurrences of the macro name itself throughout
    // the macro definition
    void replace_parameters()
    {
        using namespace boost::wave;
        
        if (!replaced_parameters) {
        typename definition_container_t::iterator end = macrodefinition.end();
        typename definition_container_t::iterator it = macrodefinition.begin(); 

            for (/**/; it != end; ++it) {
            token_id id = *it;
            
                if (T_IDENTIFIER == id || 
                    IS_CATEGORY(id, KeywordTokenType) ||
                    IS_EXTCATEGORY(id, OperatorTokenType|AltExtTokenType)) 
                {
                // may be a parameter to replace
                    const_parameter_iterator_t cend = macroparameters.end();
                    const_parameter_iterator_t cit = macroparameters.begin();
                    for (typename parameter_container_t::size_type i = 0; 
                        cit != cend; ++cit, ++i) 
                    {
                        if ((*it).get_value() == (*cit).get_value()) {
                            (*it).set_token_id(token_id(T_PARAMETERBASE+i));
                            break;
                        }
#if BOOST_WAVE_SUPPORT_VARIADICS_PLACEMARKERS != 0
                        else if (T_ELLIPSIS == token_id(*cit) && 
                            "__VA_ARGS__" == (*it).get_value()) 
                        {
                        // __VA_ARGS__ requires special handling
                            (*it).set_token_id(token_id(T_EXTPARAMETERBASE+i));
                            break;
                        }
#endif
                    }
                }
            }
            
#if BOOST_WAVE_SUPPORT_VARIADICS_PLACEMARKERS != 0 
        // we need to know, if the last of the formal arguments is an ellipsis
            if (macroparameters.size() > 0 &&
                T_ELLIPSIS == token_id(macroparameters.back())) 
            {
                has_ellipsis = true;
            }
#endif 
            replaced_parameters = true;     // do it only once
        }
    }

    TokenT macroname;                       // macro name
    parameter_container_t macroparameters;  // formal parameters
    definition_container_t macrodefinition; // macro definition token sequence
    long uid;                               // unique id of this macro
    bool is_functionlike;
    bool replaced_parameters;
    bool is_available_for_replacement;
    bool is_predefined;
#if BOOST_WAVE_SUPPORT_VARIADICS_PLACEMARKERS != 0
    bool has_ellipsis;
#endif
};

///////////////////////////////////////////////////////////////////////////////
}   // namespace util
}   // namespace wave
}   // namespace boost

#endif // !defined(MACRO_DEFINITION_HPP_D68A639E_2DA5_4E9C_8ACD_CFE6B903831E_INCLUDED)
