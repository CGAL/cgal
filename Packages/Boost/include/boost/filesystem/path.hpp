//  boost/filesystem/path.hpp  -----------------------------------------------//

//  © Copyright Beman Dawes 2002-2003
//  Use, modification, and distribution is subject to the Boost Software
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

//  See library home page at http://www.boost.org/libs/filesystem

//----------------------------------------------------------------------------// 

#ifndef BOOST_FILESYSTEM_PATH_HPP
#define BOOST_FILESYSTEM_PATH_HPP

#include <boost/filesystem/config.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <string>
#include <cassert>

#include <boost/config/abi_prefix.hpp> // must be the last header

//----------------------------------------------------------------------------//

namespace boost
{
  namespace filesystem
  {
    class directory_iterator;


  //  path -------------------------------------------------------------------//

    class BOOST_FILESYSTEM_DECL path
    {
    public:
      typedef bool (*name_check)( const std::string & name );

      // compiler generates copy constructor, copy assignment, and destructor

      path(){}

      path( const std::string & src );
      path( const char * src );

      path( const std::string & src, name_check checker );
      path( const char * src, name_check checker );

      // append operations:
      path & operator /=( const path & rhs );
      path operator /( const path & rhs ) const
        { return path( *this ) /= rhs; }

      // modification functions:
      path & normalize();

      // conversion functions:
      const std::string & string() const { return m_path; }
      std::string native_file_string() const;
      std::string native_directory_string() const;

      // decomposition functions:
      path         root_path() const;
      std::string  root_name() const;
      std::string  root_directory() const;
      path         relative_path() const;
      std::string  leaf() const;
      path         branch_path() const;

      // query functions:
      bool empty() const { return m_path.empty(); } // name consistent with std containers

      bool is_complete() const;

      bool has_root_path() const;
      bool has_root_name() const;
      bool has_root_directory() const;
      bool has_relative_path() const;
      bool has_leaf() const { return !m_path.empty(); }
      bool has_branch_path() const;

      // iteration over the names in the path:
      class iterator : public boost::iterator_facade<
        iterator,
        std::string const,
        boost::bidirectional_traversal_tag >
      {
      private:
        friend class boost::iterator_core_access;
        friend class boost::filesystem::path;

        reference dereference() const { return m_name; }
        bool equal( const iterator & rhs ) const
          { return m_path_ptr == rhs.m_path_ptr && m_pos == rhs.m_pos; }
        BOOST_FILESYSTEM_DECL void increment();
        BOOST_FILESYSTEM_DECL void decrement();

        std::string             m_name;     // cache current element.
        const path *            m_path_ptr; // path being iterated over.
        std::string::size_type  m_pos;      // position of name in
                                            // path_ptr->string(). The
                                            // end() iterator is indicated by 
                                            // pos == path_ptr->string().size()
      };

      iterator begin() const;
      iterator end() const
      {
        iterator itr;
        itr.m_path_ptr = this;
        itr.m_pos = m_path.size();
        return itr;
      }

      // default name_check mechanism:
      static bool default_name_check_writable(); 
      static void default_name_check( name_check new_check );
      static name_check default_name_check();

      // relational operators
      bool operator<( const path & that ) const;
      bool operator==( const path & that ) const { return !(*this < that) && !(that < *this); }
      bool operator!=( const path & that ) const { return !(*this == that); }
      bool operator>( const path & that ) const  { return that < *this; }
      bool operator<=( const path & that ) const { return !(that < *this); }
      bool operator>=( const path & that ) const { return !(*this < that); }

    private:
      // Note: This is an implementation for POSIX and Windows, where there
      // are only minor differences between generic and system-specific
      // constructor input formats.  Private members might be quite different
      // in other implementations, particularly where there were wide
      // differences between generic and system-specific argument formats,
      // or between native_file_string() and native_directory_string() formats.

      std::string  m_path;

      friend class directory_iterator;
      // Was qualified; como433beta8 reports:
      //    warning #427-D: qualified name is not allowed in member declaration 
      friend class iterator; 

    public: // should be private, but friend functions don't work for me
      void m_path_append( const std::string & src, name_check checker );
      void m_replace_leaf( const char * new_leaf );
    };

  //  path non-member functions  ---------------------------------------------//

    inline path operator / ( const char * lhs, const path & rhs )
      { return path( lhs ) /= rhs; }

    inline path operator / ( const std::string & lhs, const path & rhs )
      { return path( lhs ) /= rhs; }
   
  //  path::name_checks  ---------------------------------------------------//

    BOOST_FILESYSTEM_DECL bool portable_posix_name( const std::string & name );
    BOOST_FILESYSTEM_DECL bool windows_name( const std::string & name );
    BOOST_FILESYSTEM_DECL bool portable_name( const std::string & name );
    BOOST_FILESYSTEM_DECL bool portable_directory_name( const std::string & name );
    BOOST_FILESYSTEM_DECL bool portable_file_name( const std::string & name );
    BOOST_FILESYSTEM_DECL bool no_check( const std::string & name );   // always returns true
    BOOST_FILESYSTEM_DECL bool native( const std::string & name );
      // native(name) must return true for any name which MIGHT be valid
      // on the native platform.

  } // namespace filesystem
} // namespace boost

#include <boost/config/abi_suffix.hpp> // pops abi_suffix.hpp pragmas
#endif // BOOST_FILESYSTEM_PATH_HPP
