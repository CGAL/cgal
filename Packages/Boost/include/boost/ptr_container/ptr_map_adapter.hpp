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

#ifndef BOOST_PTR_CONTAINER_DETAIL_PTR_MAP_ADAPTER_HPP
#define BOOST_PTR_CONTAINER_DETAIL_PTR_MAP_ADAPTER_HPP

#if defined(_MSC_VER) && (_MSC_VER >= 1200)
# pragma once
#endif

#include <boost/ptr_container/detail/map_iterator.hpp>
#include <boost/ptr_container/detail/associative_ptr_container.hpp>
#include <boost/static_assert.hpp>
#include <boost/range/iterator_range.hpp>

namespace boost
{
namespace ptr_container_detail
{

    template
    < 
        class T,
        class VoidPtrMap
    >
    struct map_config
    {
        typedef BOOST_DEDUCED_TYPENAME remove_nullable<T>::type
                     U;
        typedef VoidPtrMap 
                     void_container_type;
        
        typedef BOOST_DEDUCED_TYPENAME VoidPtrMap::allocator_type
                     allocator_type;
        
        typedef BOOST_DEDUCED_TYPENAME VoidPtrMap::key_compare
                     key_compare;
        
        typedef BOOST_DEDUCED_TYPENAME VoidPtrMap::value_compare
                     value_compare;
        
        typedef BOOST_DEDUCED_TYPENAME VoidPtrMap::key_type
                     key_type;
        
        typedef U    value_type;

        typedef ptr_map_iterator< 
                       BOOST_DEDUCED_TYPENAME VoidPtrMap::iterator,
                       BOOST_DEDUCED_TYPENAME VoidPtrMap::key_type, value_type>
                     iterator;
        
        typedef
            ptr_map_iterator<
                       BOOST_DEDUCED_TYPENAME VoidPtrMap::const_iterator,
                       BOOST_DEDUCED_TYPENAME VoidPtrMap::key_type, 
                       const value_type>
                     const_iterator;

        typedef
            ptr_map_iterator<
                       BOOST_DEDUCED_TYPENAME VoidPtrMap::reverse_iterator,
                       BOOST_DEDUCED_TYPENAME VoidPtrMap::key_type, value_type>
                     reverse_iterator;

        typedef
            ptr_map_iterator<
                       BOOST_DEDUCED_TYPENAME VoidPtrMap::const_reverse_iterator,
                       BOOST_DEDUCED_TYPENAME VoidPtrMap::key_type, 
                       const value_type>
                     const_reverse_iterator;

        typedef std::pair<const key_type, void*>
                     object_type;

        template< class Iter >
        static U* get_pointer( Iter i )
        {
            return static_cast<U*>( i.base()->second );
        }

        template< class Iter >
        static const U* get_const_pointer( Iter i )
        {
            return static_cast<const U*>( i.base()->second );
        }

        BOOST_STATIC_CONSTANT( bool, allow_null = boost::is_nullable<T>::value );
    };
    
    

    template
    < 
        class T,
        class VoidPtrMap, 
        class CloneAllocator
    >
    class ptr_map_adapter_base : 
        public ptr_container_detail::associative_ptr_container< map_config<T,VoidPtrMap>,
                                                    CloneAllocator >
    {
        typedef ptr_container_detail::associative_ptr_container< map_config<T,VoidPtrMap>,
                                                     CloneAllocator > 
            base_type;

        typedef ptr_map_adapter_base<T,VoidPtrMap,CloneAllocator> this_type;
        
    public:

        typedef BOOST_DEDUCED_TYPENAME base_type::const_iterator
                    const_iterator;
        typedef BOOST_DEDUCED_TYPENAME base_type::key_type
                    key_type;
        typedef BOOST_DEDUCED_TYPENAME base_type::reference
                    reference;
        typedef BOOST_DEDUCED_TYPENAME base_type::value_type
                    value_type;
        typedef BOOST_DEDUCED_TYPENAME base_type::auto_type
                    auto_type;

    private:
        reference lookup( const key_type& key ) const
        {
           iterator i = const_cast<ptr_map_adapter_base*>(this)
                                          ->find( key );
           if( i != const_cast<ptr_map_adapter_base*>(this)->end() )
               return *i;
           else                                           
               throw bad_ptr_container_operation( "'ptr_map/multimap::at()' could"
                                                    " not find key" );
        }

        struct eraser // scope guard
        {
            bool            released_;
            VoidPtrMap*     m_;
            const key_type& key_;

            eraser( VoidPtrMap* m, const key_type& key ) 
              : released_(false), m_(m), key_(key)
            {}

            ~eraser() 
            {
                if( !released_ )
                    m_->erase(key_);
            }

            void release() { released_ = true; }
        };

        reference insert_lookup( const key_type& key )
        {
            void*& ref = this->c_private()[key];
            if( ref )
            {
                value_type v = static_cast<value_type>( ref );
                return *v;
            }
            else
            {
                eraser e(&this->c_private(),key); // nothrow
                value_type res = new T();         // strong 
                ref = res;                        // nothrow
                e.release();                      // nothrow
                return *res;
            }
          }
        
    public:

        BOOST_PTR_CONTAINER_DEFINE_CONSTRUCTORS( ptr_map_adapter_base, 
                                                 base_type );

        template< class Compare, class Allocator >
        explicit ptr_map_adapter_base( const Compare& comp,
                                       const Allocator& a ) 
        : base_type( comp, a ) 
        { }

        template< class PtrContainer >
        ptr_map_adapter_base( std::auto_ptr<PtrContainer> clone ) 
        : base_type( clone )
        { }
        
        template< typename PtrContainer >
        void operator=( std::auto_ptr<PtrContainer> clone )    
        {
            base_type::operator=( clone );
        }        

        iterator find( const key_type& x )                                                
        {                                                                            
            return iterator( this->c_private().find( x ) );                                
        }                                                                            

        const_iterator find( const key_type& x ) const                                    
        {                                                                            
            return const_iterator( this->c_private().find( x ) );                          
        }                                                                            

        size_type count( const key_type& x ) const                                        
        {                                                                            
            return this->c_private().count( x );                                           
        }                                                                            
                                                                                     
        iterator lower_bound( const key_type& x )                                         
        {                                                                            
            return iterator( this->c_private().lower_bound( x ) );                         
        }                                                                            
                                                                                     
        const_iterator lower_bound( const key_type& x ) const                             
        {                                                                            
            return const_iterator( this->c_private().lower_bound( x ) );                   
        }                                                                            
                                                                                     
        iterator upper_bound( const key_type& x )                                         
        {                                                                            
            return iterator( this->c_private().upper_bound( x ) );                         
        }                                                                            
                                                                                     
        const_iterator upper_bound( const key_type& x ) const                             
        {                                                                            
            return const_iterator( this->c_private().upper_bound( x ) );                   
        }                                                                            
                                                                                     
        iterator_range<iterator> equal_range( const key_type& x )                    
        {                                                                            
            std::pair<BOOST_DEDUCED_TYPENAME base_type::ptr_iterator,
                      BOOST_DEDUCED_TYPENAME base_type::ptr_iterator>
                 p = this->c_private().equal_range( x );   
            return make_iterator_range( iterator( p.first ), iterator( p.second ) );      
        }                                                                            
                                                                                     
        iterator_range<const_iterator> equal_range( const key_type& x ) const  
        {                                                                            
            std::pair<BOOST_DEDUCED_TYPENAME base_type::ptr_const_iterator,
                      BOOST_DEDUCED_TYPENAME base_type::ptr_const_iterator> 
                p = this->c_private().equal_range( x ); 
            return make_iterator_range( const_iterator( p.first ), const_iterator( p.second ) );    
        }                                                                            
                                                                                     
        reference at( const key_type& key )                                              
        {                   
            return lookup( key );                                                         
        }                                                                            
                                                                                     
        const_reference at( const key_type& key ) const                                  
        {                                                                            
            return lookup( key );
        }

        reference operator[]( const key_type& key )                                              
        {                          
            return insert_lookup( key );
        }              

        auto_type replace( iterator where, value_type x ) // strong  
        { 
            BOOST_ASSERT( where != this->end() );

            this->enforce_null_policy( x, "Null pointer in 'replace()'" );

            auto_type ptr( x );

            if( this->empty() )
                throw bad_ptr_container_operation( "'replace()' on empty container" );

            auto_type old( &*where );               // nothrow
            where.base()->second = ptr.release();   // nothrow, commit
            return move( old );
        }
                                                                                     
    };
    
} // ptr_container_detail

    /////////////////////////////////////////////////////////////////////////
    // ptr_map_adapter
    /////////////////////////////////////////////////////////////////////////
    
    template
    < 
        class T,
        class VoidPtrMap, 
        class CloneAllocator = heap_clone_allocator
    >
    class ptr_map_adapter : 
        public ptr_container_detail::ptr_map_adapter_base<T,VoidPtrMap,CloneAllocator>
    {
        typedef ptr_container_detail::ptr_map_adapter_base<T,VoidPtrMap,CloneAllocator> 
            base_type;
    
    public:    
        typedef BOOST_DEDUCED_TYPENAME base_type::iterator 
                     iterator;       
        typedef BOOST_DEDUCED_TYPENAME base_type::const_iterator
                     const_iterator;
        typedef BOOST_DEDUCED_TYPENAME base_type::object_type
                    object_type;
        typedef BOOST_DEDUCED_TYPENAME base_type::size_type
                    size_type;
        typedef BOOST_DEDUCED_TYPENAME base_type::key_type
                    key_type;
        typedef BOOST_DEDUCED_TYPENAME base_type::const_reference
                    const_reference;
        typedef BOOST_DEDUCED_TYPENAME base_type::auto_type
                    auto_type;
        typedef BOOST_DEDUCED_TYPENAME VoidPtrMap::key_compare 
                    key_compare;
        typedef BOOST_DEDUCED_TYPENAME VoidPtrMap::allocator_type 
                    allocator_type;
        typedef BOOST_DEDUCED_TYPENAME base_type::value_type
                    value_type;
    private:

        void safe_insert( const key_type& key, auto_type ptr ) // strong
        {
            std::pair<BOOST_DEDUCED_TYPENAME base_type::ptr_iterator,bool>
                res = 
                this->c_private().insert( std::make_pair( key, ptr.get() ) ); // strong, commit      
            if( res.second )                                                  // nothrow     
                ptr.release();                                                // nothrow
        }
    
        template< class II >                                               
        void map_basic_clone_and_insert( II first, II last )                  
        {       
            while( first != last )                                            
            {                                            
                if( this->find( first.key() ) == this->end() )
                {
                    const object_type& p = *first.base();     // nothrow                    
                    auto_type ptr( this->null_policy_allocate_clone(
                                static_cast<value_type>(p.second) ) ); 
                                                              // strong 
                    this->safe_insert( p.first, ptr_container_detail::
                                                move( ptr ) );// strong, commit
                }
                ++first;                                                      
            }                                                                 
        }
    
    public:

        explicit ptr_map_adapter( const key_compare& comp = key_compare(),
                                  const allocator_type& a = allocator_type()  ) 
          : base_type( comp, a ) { }
    
        template< class InputIterator >
        ptr_map_adapter( InputIterator first, InputIterator last, 
                         const key_compare& comp = key_compare(),
                         const allocator_type& a = allocator_type() )
          : base_type( comp, a ) 
        {
            map_basic_clone_and_insert( first, last );
        }

        template< class U >
        ptr_map_adapter( std::auto_ptr<U> r ) : base_type( r )
        { }

        template< class U >
        void operator=( std::auto_ptr<U> r )
        {  
            base_type::operator=( r );
        }

        using base_type::release;

        template< typename InputIterator >
        void insert( InputIterator first, InputIterator last ) // basic
        {
            map_basic_clone_and_insert( first, last );
        }

        template< class Range >
        void insert( const Range& r )
        {
            insert( this->adl_begin(r), this->adl_end(r) );
        }

        std::pair<iterator,bool> insert( key_type& key, value_type x ) // strong
        {
            this->enforce_null_policy( x, "Null pointer in ptr_map_adapter::insert()" );
            auto_type ptr( x );                                                 // nothrow
    
            std::pair<BOOST_DEDUCED_TYPENAME base_type::ptr_iterator,bool>
                 res = this->c_private().insert( std::make_pair( key, x ) );       // strong, commit      
            if( res.second )                                                                               // nothrow     
                ptr.release();                                                                             // nothrow
            return std::make_pair( iterator( res.first ), res.second );  // nothrow   
        }

        bool transfer( iterator object, 
                       ptr_map_adapter& from ) // strong
        {
            return this->single_transfer( object, from );
        }

        size_type transfer( iterator first, 
                            iterator last, 
                            ptr_map_adapter& from ) // basic
        {
            return this->single_transfer( first, last, from );
        }

#ifdef BOOST_NO_SFINAE
#else    

        template< class Range >
        BOOST_DEDUCED_TYPENAME boost::disable_if< boost::is_same< Range,
                                                                  iterator >,
                                                            size_type >::type
        transfer( const Range& r, ptr_map_adapter& from ) // basic
        {
            return transfer( this->adl_begin(r), this->adl_end(r), from );
        }
        
#endif

        size_type transfer( ptr_map_adapter& from ) // basic
        {
            return transfer( from.begin(), from.end(), from );
        }
        
  };
  
  /////////////////////////////////////////////////////////////////////////
  // ptr_multimap_adapter
  /////////////////////////////////////////////////////////////////////////

    template
    < 
        class T,
        class VoidPtrMultiMap, 
        class CloneAllocator = heap_clone_allocator
    >
    class ptr_multimap_adapter : 
        public ptr_container_detail::ptr_map_adapter_base<T,VoidPtrMultiMap,CloneAllocator>
    {
        typedef ptr_container_detail::ptr_map_adapter_base<T,VoidPtrMultiMap,CloneAllocator>
             base_type;
        
    public: // typedefs
        typedef BOOST_DEDUCED_TYPENAME base_type::iterator           
                       iterator;                 
        typedef BOOST_DEDUCED_TYPENAME base_type::const_iterator     
                       const_iterator;           
        typedef BOOST_DEDUCED_TYPENAME base_type::object_type
                       object_type;         
        typedef BOOST_DEDUCED_TYPENAME base_type::size_type
                       size_type;
        typedef BOOST_DEDUCED_TYPENAME base_type::key_type
                       key_type;
        typedef BOOST_DEDUCED_TYPENAME base_type::const_reference
                       const_reference;
        typedef BOOST_DEDUCED_TYPENAME base_type::value_type
                    value_type;
        typedef BOOST_DEDUCED_TYPENAME base_type::auto_type
                    auto_type;            
        typedef BOOST_DEDUCED_TYPENAME VoidPtrMultiMap::key_compare 
                    key_compare;
        typedef BOOST_DEDUCED_TYPENAME VoidPtrMultiMap::allocator_type 
                    allocator_type;
    private:

        void safe_insert( const key_type& key, auto_type ptr ) // strong
        {
            this->c_private().insert( 
                           std::make_pair( key, ptr.get() ) ); // strong, commit      
            ptr.release();                                     // nothrow
        }
        
        template< typename II >                                               
        void map_basic_clone_and_insert( II first, II last )                  
        {                                                         
            while( first != last )                                            
            {                                            
                const object_type& pair = *first.base();  // nothrow                     
                auto_type ptr( this->null_policy_allocate_clone(
                                static_cast<value_type>( pair.second ) ) );    
                                                          // strong
                safe_insert( pair.first, ptr_container_detail::
                                         move( ptr ) );   // strong, commit
                ++first;                                                      
            }                                                                 
        }
        
    public:
        
        explicit ptr_multimap_adapter( const key_compare& comp = key_compare(),
                                       const allocator_type& a = allocator_type() )
          : base_type( comp, a ) { }
        
        template< class InputIterator >
        ptr_multimap_adapter( InputIterator first, InputIterator last,
                              const key_compare& comp = key_compare(),
                              const allocator_type& a = allocator_type() )
          : base_type( comp, a )
        {
            map_basic_clone_and_insert( first, last );
        }

        template< class U >
        ptr_multimap_adapter( std::auto_ptr<U> r ) : base_type( r )
        { }

        template< class U >
        void operator=( std::auto_ptr<U> r )
        {  
            base_type::operator=( r );
        }

        using base_type::release;
        
        template< typename InputIterator >
        void insert( InputIterator first, InputIterator last ) // basic
        {
            map_basic_clone_and_insert( first, last );
        }

        template< class Range >
        void insert( const Range& r )
        {
            insert( this->adl_begin(r), this->adl_end(r) );
        }

        iterator insert( key_type& key, value_type x ) // strong
        {
            this->enforce_null_policy( x, "Null pointer in 'ptr_multimap_adapter::insert()'" );

            auto_type ptr( x );         // nothrow
            BOOST_DEDUCED_TYPENAME base_type::ptr_iterator
                res = this->c_private().insert( std::make_pair( key, x ) );
                                        // strong, commit        
            ptr.release();              // notrow
            return iterator( res );           
        }

        
        void transfer( iterator object, 
                       ptr_multimap_adapter& from ) // strong
        {
            this->multi_transfer( object, from );
        }

        size_type transfer( iterator first, 
                            iterator last, 
                            ptr_multimap_adapter& from ) // basic
        {
            return this->multi_transfer( first, last, from );
        }

#ifdef BOOST_NO_SFINAE
#else    

        template< class Range >
        BOOST_DEDUCED_TYPENAME boost::disable_if< boost::is_same< Range,
                                                                  iterator >,
                                                            size_type >::type
        transfer(  const Range& r, ptr_multimap_adapter& from ) // basic
        {
            return transfer( this->adl_begin(r), this->adl_end(r), from );
        }

#endif        

        void transfer( ptr_multimap_adapter& from ) // basic
        {
            transfer( from.begin(), from.end(), from );
            BOOST_ASSERT( from.empty() );
        }
    };

    template< class I, class K, class V >
    inline bool is_null( ptr_map_iterator<I,K,V> i )
    {
        return i.base()->second == 0;
    }
    
} // namespace 'boost'  

#endif
