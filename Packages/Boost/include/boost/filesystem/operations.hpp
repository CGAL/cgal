//  boost/filesystem/directory.hpp  ------------------------------------------//

//  Copyright © 2002, 2003 Beman Dawes
//  Copyright © 2002 Jan Langer
//  Copyright © 2001 Dietmar Kühl                                        
//  
//  Use, modification, and distribution is subject to the Boost Software
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy
//  at http://www.boost.org/LICENSE_1_0.txt)                             

//  See library home page at http://www.boost.org/libs/filesystem

//----------------------------------------------------------------------------// 

#ifndef BOOST_FILESYSTEM_DIRECTORY_HPP
#define BOOST_FILESYSTEM_DIRECTORY_HPP

#include <boost/filesystem/path.hpp>  // includes <boost/filesystem/config.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/iterator.hpp>
#include <boost/cstdint.hpp>

#include <string>
#include <ctime>

#include <boost/config/abi_prefix.hpp> // must be the last header

# ifdef BOOST_NO_STDC_NAMESPACE
    namespace std { using ::time_t; }
# endif

//----------------------------------------------------------------------------//

namespace boost
{
  namespace filesystem
  {

//  query functions  ---------------------------------------------------------//

    BOOST_FILESYSTEM_DECL bool exists( const path & ph );
    BOOST_FILESYSTEM_DECL bool symbolic_link_exists( const path & ph );
    BOOST_FILESYSTEM_DECL bool is_directory( const path & ph );

    // VC++ 7.0 and earlier has a serious namespace bug that causes a clash
    // between boost::filesystem::is_empty and the unrelated type trait
    // boost::is_empty. The workaround for those who must use broken versions
    // of VC++ is to use the function _is_empty. All others should use the
    // correct is_empty name.
    BOOST_FILESYSTEM_DECL bool _is_empty( const path & ph ); // deprecated

#   if !defined( BOOST_MSVC ) || BOOST_MSVC > 1300
    inline bool is_empty( const path & ph ) { return _is_empty( ph ); }
#   endif

    BOOST_FILESYSTEM_DECL bool equivalent( const path & ph1, const path & ph2 );
    BOOST_FILESYSTEM_DECL boost::intmax_t file_size( const path & ph );
    BOOST_FILESYSTEM_DECL std::time_t last_write_time( const path & ph );
    BOOST_FILESYSTEM_DECL void last_write_time( const path & ph, const std::time_t new_time );

//  operations  --------------------------------------------------------------//

    BOOST_FILESYSTEM_DECL bool create_directory( const path & directory_ph );

    BOOST_FILESYSTEM_DECL bool remove( const path & ph );
    BOOST_FILESYSTEM_DECL unsigned long remove_all( const path & ph );

    BOOST_FILESYSTEM_DECL void rename( const path & from_path,
                 const path & to_path );

    BOOST_FILESYSTEM_DECL void copy_file( const path & from_file_ph,
                    const path & to_file_ph );

    BOOST_FILESYSTEM_DECL path current_path();
    BOOST_FILESYSTEM_DECL const path & initial_path();

    BOOST_FILESYSTEM_DECL path system_complete( const path & ph );
    BOOST_FILESYSTEM_DECL path complete( const path & ph, const path & base = initial_path() );

//  test helper  -------------------------------------------------------------//

    // not part of the documented interface because false positives are possible;
    // there is no law that says that an OS that has large stat.st_size
    // actually supports large file sizes.
    BOOST_FILESYSTEM_DECL bool possible_large_file_size_support();


//  directory_iterator helpers  ----------------------------------------------//
//    forwarding functions avoid need for BOOST_FILESYSTEM_DECL for class
//    directory_iterator, and so avoid iterator_facade DLL template problems
    namespace detail
    {
      class dir_itr_imp;
      // shared_ptr provides shallow-copy semantics required for InputIterators
      typedef boost::shared_ptr< dir_itr_imp > dir_itr_imp_ptr;
      BOOST_FILESYSTEM_DECL void dir_itr_init( dir_itr_imp_ptr & m_imp,
                                               const path & dir_path );
      BOOST_FILESYSTEM_DECL path & dir_itr_dereference(
                                               const dir_itr_imp_ptr & m_imp );
      BOOST_FILESYSTEM_DECL void dir_itr_increment( dir_itr_imp_ptr & m_imp );
    } // namespace detail

//  directory_iterator  ------------------------------------------------------//

    class directory_iterator
      : public boost::iterator_facade<
      directory_iterator,
      path,
      boost::single_pass_traversal_tag >
    {
    public:
      directory_iterator(){}  // creates the "end" iterator
      explicit directory_iterator( const path & p )
        { detail::dir_itr_init( m_imp, p ); }

/*
The *r++ requirement doesn't appear to apply to the new single_pass_traversal category
Thus I'm leaving the proxy out pending confirmation from the N1477 authors
struct path_proxy // allows *r++ to work, as required by 24.1.1
      {
        path pv;
        explicit path_proxy( const path & p ) : pv(p) {}
        path operator*() const { return pv; }
      };

      path_proxy operator++(int)
      {
        path_proxy pp( m_deref() );
        ++*this;
        return pp;
      }
*/

    private:
      detail::dir_itr_imp_ptr  m_imp;
      friend class boost::iterator_core_access;
      reference dereference() const 
        { return detail::dir_itr_dereference( m_imp ); }
      void increment()
        { detail::dir_itr_increment( m_imp ); }
      bool equal( const directory_iterator & rhs ) const
        { return m_imp == rhs.m_imp; }
    };
  } // namespace filesystem
} // namespace boost


#include <boost/config/abi_suffix.hpp> // pops abi_suffix.hpp pragmas
#endif // BOOST_FILESYSTEM_DIRECTORY_HPP
