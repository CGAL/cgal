//
// Boost.Pointer Container
//
//  Copyright Thorsten Ottosen 2003-2005. Use, modification and
//  distribution is subject to the Boost Software License, Version
//  1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
// For more information, see http://www.boost.org/libs/ptr_container/
//

#ifndef BOOST_PTR_CONTAINER_CLONE_ALLOCATOR_HPP
#define BOOST_PTR_CONTAINER_CLONE_ALLOCATOR_HPP

#include <boost/checked_delete.hpp>

namespace boost
{
    /////////////////////////////////////////////////////////////////////////
    // Clonable concept 
    /////////////////////////////////////////////////////////////////////////
    
    template< class T >
    inline T* new_clone( const T& r )
    {
        return new T( r );
    }

    template< class T >
    inline void delete_clone( const T* r )
    {
        checked_delete( r );
    }

    /////////////////////////////////////////////////////////////////////////
    // CloneAllocator concept
    /////////////////////////////////////////////////////////////////////////
    
    struct heap_clone_allocator
    {
        template< class U >
        static U* allocate_clone( const U& r )
        {
            return new_clone( r );
        }

        template< class U >
        static void deallocate_clone( const U* r )
        {
            delete_clone( r );
        }

    };


    
    struct view_clone_allocator
    {
        template< class U >
        static U* allocate_clone( const U& r )
        {
            return const_cast<U*>( &r );
        }

        template< class U >
        static void deallocate_clone( const U* r )
        {
            // do nothing
        }
    };

    /////////////////////////////////////////////////////////////////////////
    // MapCloneAllocator concept
    /////////////////////////////////////////////////////////////////////////

    template< class T >
    inline T* new_default_clone( const T* )
    {
        return new T();
    }
        
    struct map_heap_clone_allocator : heap_clone_allocator
    {   
        template< class U >
        static U* allocate_default_clone()
        {
            static const U* ptr = 0;
            return new_default_clone(ptr);
        }
    };
    
} // namespace 'boost'

#endif

