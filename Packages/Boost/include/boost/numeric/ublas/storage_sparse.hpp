//
//  Copyright (c) 2000-2002
//  Joerg Walter, Mathias Koch
//
//  Permission to use, copy, modify, distribute and sell this software
//  and its documentation for any purpose is hereby granted without fee,
//  provided that the above copyright notice appear in all copies and
//  that both that copyright notice and this permission notice appear
//  in supporting documentation.  The authors make no representations
//  about the suitability of this software for any purpose.
//  It is provided "as is" without express or implied warranty.
//
//  The authors gratefully acknowledge the support of
//  GeNeSys mbH & Co. KG in producing this work.
//

#ifndef BOOST_UBLAS_STORAGE_SPARSE_H
#define BOOST_UBLAS_STORAGE_SPARSE_H

#include <map>

#include <boost/numeric/ublas/storage.hpp>

namespace boost { namespace numeric { namespace ublas {

    namespace detail {

        template<class I, class T, class C>
        BOOST_UBLAS_INLINE
        I lower_bound (const I &begin, const I &end, const T &t, C compare) {
            // t <= *begin <=> ! (*begin < t)
            if (begin == end || ! compare (*begin, t))
                return begin;
            if (compare (*(end - 1), t))
                return end;
            return std::lower_bound (begin, end, t, compare);
        }
        template<class I, class T, class C>
        BOOST_UBLAS_INLINE
        I upper_bound (const I &begin, const I &end, const T &t, C compare) {
            if (begin == end || compare (t, *begin))
                return begin;
            // (*end - 1) <= t <=> ! (t < *end)
            if (! compare (t, *(end - 1)))
                return end;
            return std::upper_bound (begin, end, t, compare);
        }

        template<class P>
        struct less_pair {
            BOOST_UBLAS_INLINE
            bool operator () (const P &p1, const P &p2) {
                return p1.first < p2.first;
            }
        };
        template<class T>
        struct less_triple {
            BOOST_UBLAS_INLINE
            bool operator () (const T &t1, const T &t2) {
                return t1.first.first < t2.first.first ||
                       (t1.first.first == t2.first.first && t1.first.second < t2.first.second);
            }
        };

    }

#ifdef BOOST_UBLAS_STRICT_MAP_ARRAY
    template<class A>
    class sparse_storage_element:
       public container_reference<A> {
    public:
        typedef A array_type;
        typedef typename A::key_type index_type;
        typedef typename A::mapped_type data_value_type;
        // typedef const data_value_type &data_const_reference;
        typedef typename type_traits<data_value_type>::const_reference data_const_reference;
        typedef data_value_type &data_reference;
        typedef typename A::value_type value_type;
        typedef value_type *pointer;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        sparse_storage_element (array_type &a, pointer it):
            container_reference<array_type> (a), it_ (it), i_ (it->first), d_ (it->second), dirty_ (false) {}
        BOOST_UBLAS_INLINE
        sparse_storage_element (array_type &a, index_type i):
            container_reference<array_type> (a), it_ (), i_ (i), d_ (), dirty_ (false) {
            pointer it = (*this) ().find (i_);
            if (it == (*this) ().end ())
                it = (*this) ().insert ((*this) ().end (), value_type (i_, d_));
            d_ = it->second;
        }
        BOOST_UBLAS_INLINE
        ~sparse_storage_element () {
            if (dirty_) {
                if (! it_)
                    it_ = (*this) ().find (i_);
                BOOST_UBLAS_CHECK (it_ != (*this) ().end (), internal_logic ());
                it_->second = d_;
            }
        }

        // Element access - only if data_const_reference is defined
        BOOST_UBLAS_INLINE
        typename data_value_type::data_const_reference
        operator [] (index_type i) const {
            return d_ [i];
        }

        // Assignment
        BOOST_UBLAS_INLINE
        sparse_storage_element &operator = (const sparse_storage_element &p) {
            // Overide the implict copy assignment
            d_ = p.d_;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_storage_element &operator = (const D &d) {
            d_ = d;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_storage_element &operator += (const D &d) {
            d_ += d;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_storage_element &operator -= (const D &d) {
            d_ -= d;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_storage_element &operator *= (const D &d) {
            d_ *= d;
            dirty_ = true;
            return *this;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        sparse_storage_element &operator /= (const D &d) {
            d_ /= d;
            dirty_ = true;
            return *this;
        }

        // Comparison
        template<class D>
        BOOST_UBLAS_INLINE
        bool operator == (const D &d) const {
            return d_ == d;
        }
        template<class D>
        BOOST_UBLAS_INLINE
        bool operator != (const D &d) const {
            return d_ != d;
        }

        // Conversion
        BOOST_UBLAS_INLINE
        operator data_const_reference () const {
            return d_;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (sparse_storage_element p) {
            if (this != &p) {
                dirty_ = true;
                p.dirty_ = true;
                std::swap (d_, p.d_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (sparse_storage_element p1, sparse_storage_element p2) {
            p1.swap (p2);
        }
#endif

    private:
        pointer it_;
        index_type i_;
        data_value_type d_;
        bool dirty_;
    };
#endif


    // Default map type is simply forwarded to std::map
    // FIXME should use ALLOC for map but std::pair<const I, T> and std::pair<I,T> fail
    template<class I, class T, class ALLOC>
    class map_std : public std::map<I, T /*, ALLOC */> {
    };


    // Map array
    //  Implementation requires pair<I, T> allocator definition (without const)
    template<class I, class T, class ALLOC>
    class map_array {
    public:
        typedef ALLOC allocator_type;
        typedef typename ALLOC::size_type size_type;
        typedef typename ALLOC::difference_type difference_type;
        typedef std::pair<I,T> value_type;
        typedef I key_type;
        typedef T mapped_type;
        typedef const value_type &const_reference;
        typedef value_type &reference;
        typedef const value_type *const_pointer;
        typedef value_type *pointer;
        typedef const T &data_const_reference;
#ifndef BOOST_UBLAS_STRICT_MAP_ARRAY
        typedef T &data_reference;
#else
        typedef sparse_storage_element<map_array> data_reference;
#endif

        // Construction and destruction
        BOOST_UBLAS_INLINE
        map_array (const ALLOC &a = ALLOC()):
            alloc_(a), capacity_ (0), size_ (0) {
                data_ = 0;
        }
        BOOST_UBLAS_INLINE
        map_array (const map_array &c):
            alloc_ (c.alloc_), capacity_ (c.size_), size_ (c.size_) {
            if (capacity_) {
                data_ = alloc_.allocate (capacity_ BOOST_UBLAS_ALLOCATOR_HINT);
                std::uninitialized_copy (data_, data_ + capacity_, c.data_);
                // capacity != size_ requires uninitialized_fill (size_ to capacity_)
            }
            else
                data_ = 0;
        }
        BOOST_UBLAS_INLINE
        ~map_array () {
            if (capacity_) {
                std::for_each (data_, data_ + capacity_, static_destroy);
                alloc_.deallocate (data_, capacity_);
            }
        }

    private:
        // Resizing - implicitly exposses uninitialized (but default constructed) mapped_type
        BOOST_UBLAS_INLINE
        void resize (size_type size) {
            BOOST_UBLAS_CHECK (size_ <= capacity_, internal_logic ());
            if (size > capacity_) {
                const size_type capacity = size << 1;
                BOOST_UBLAS_CHECK (capacity, internal_logic ());
                pointer data = alloc_.allocate (capacity BOOST_UBLAS_ALLOCATOR_HINT);
                std::uninitialized_copy (data_, data_ + (std::min) (size, size_), data);
                std::uninitialized_fill (data + (std::min) (size, size_), data + capacity, value_type ());

                if (capacity_) {
                    std::for_each (data_, data_ + capacity_, static_destroy);
                    alloc_.deallocate (data_, capacity_);
                }
                capacity_ = capacity;
                data_ = data;
            }
            size_ = size;
            BOOST_UBLAS_CHECK (size_ <= capacity_, internal_logic ());
        }
    public:

        // Reserving
        BOOST_UBLAS_INLINE
        void reserve (size_type capacity) {
            BOOST_UBLAS_CHECK (size_ <= capacity_, internal_logic ());
            // Reduce capacity_ if size_ allows
            BOOST_UBLAS_CHECK (capacity >= size_, bad_size ());
            pointer data;
            if (capacity) {
                data = alloc_.allocate (capacity BOOST_UBLAS_ALLOCATOR_HINT);
                std::uninitialized_copy (data_, data_ + size_, data);
                std::uninitialized_fill (data + size_, data + capacity, value_type ());
            }
            else
                data = 0;
                
            if (capacity_) {
                std::for_each (data_, data_ + capacity_, static_destroy);
                alloc_.deallocate (data_, capacity_);
            }
            capacity_ = capacity;
            data_ = data;
            BOOST_UBLAS_CHECK (size_ <= capacity_, internal_logic ());
        }

        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }
        BOOST_UBLAS_INLINE
        size_type capacity () const {
            return capacity_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        data_reference operator [] (key_type i) {
#ifndef BOOST_UBLAS_STRICT_MAP_ARRAY
            pointer it = find (i);
            if (it == end ())
                it = insert (end (), value_type (i, mapped_type (0)));
            BOOST_UBLAS_CHECK (it != end (), internal_logic ());
            return it->second;
#else
            return data_reference (*this, i);
#endif
        }

        // Assignment
        BOOST_UBLAS_INLINE
        map_array &operator = (const map_array &a) {
            if (this != &a) {
                resize (a.size_);
                std::copy (a.data_, a.data_ + a.size_, data_);
            }
            return *this;
        }
        BOOST_UBLAS_INLINE
        map_array &assign_temporary (map_array &a) {
            swap (a);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (map_array &a) {
            if (this != &a) {
                std::swap (capacity_, a.capacity_);
                std::swap (data_, a.data_);
                std::swap (size_, a.size_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (map_array &a1, map_array &a2) {
            a1.swap (a2);
        }
#endif

        // Element insertion and deletion
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer push_back (pointer it, const value_type &p) {
            if (size () == 0 || (it = end () - 1)->first < p.first) {
                resize (size () + 1);
                *(it = end () - 1) = p;
                return it;
            }
            external_logic ().raise ();
            return it;
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer insert (pointer it, const value_type &p) {
            it = detail::lower_bound (begin (), end (), p, detail::less_pair<value_type> ());
            difference_type n = it - begin ();
            BOOST_UBLAS_CHECK (size () == 0 || size () == size_type (n) || it->first != p.first, external_logic ());
            resize (size () + 1);
            it = begin () + n;
            std::copy_backward (it, end () - 1, end ());
            *it = p;
            return it;
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        void erase (pointer it) {
            BOOST_UBLAS_CHECK (begin () <= it && it < end (), bad_index ());
            // Fixed by George Katsirelos.
            // (*it).second = mapped_type (0);
            std::copy (it + 1, end (), it);
            resize (size () - 1);
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        void erase (pointer it1, pointer it2) {
            BOOST_UBLAS_CHECK (begin () <= it1 && it1 < it2 && it2 <= end (), bad_index ());
            // Fixed by George Katsirelos.
            // while (it1 != it2) {
            //     BOOST_UBLAS_CHECK (begin () <= it1 && it1 < end (), bad_index ());
            //     (*it1).second = mapped_type (0);
            //     ++ it1;
            // }
            std::copy (it2, end (), it1);
            resize (size () - (it2 - it1));
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        void clear () {
            resize (0);
        }

        // Element lookup
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_pointer find (key_type i) const {
            const_pointer it (detail::lower_bound (begin (), end (), value_type (i, mapped_type (0)), detail::less_pair<value_type> ()));
            if (it == end () || it->first != i)
                it = end ();
            return it;
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer find (key_type i) {
            pointer it (detail::lower_bound (begin (), end (), value_type (i, mapped_type (0)), detail::less_pair<value_type> ()));
            if (it == end () || it->first != i)
                it = end ();
            return it;
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_pointer lower_bound (key_type i) const {
            return detail::lower_bound (begin (), end (), value_type (i, mapped_type (0)), detail::less_pair<value_type> ());
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer lower_bound (key_type i) {
            return detail::lower_bound (begin (), end (), value_type (i, mapped_type (0)), detail::less_pair<value_type> ());
        }

        // Iterators simply are pointers.

        typedef const_pointer const_iterator;

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return data_;
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return data_ + size_;
        }

        typedef pointer iterator;

        BOOST_UBLAS_INLINE
        iterator begin () {
            return data_;
        }
        BOOST_UBLAS_INLINE
        iterator end () {
            return data_ + size_;
        }

        // Reverse iterators

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef std::reverse_iterator<const_iterator, value_type, const_reference> const_reverse_iterator;
#else
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
#endif

        BOOST_UBLAS_INLINE
        const_reverse_iterator rbegin () const {
            return const_reverse_iterator (end ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator rend () const {
            return const_reverse_iterator (begin ());
        }

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef std::reverse_iterator<iterator, value_type, reference> reverse_iterator;
#else
        typedef std::reverse_iterator<iterator> reverse_iterator;
#endif

        BOOST_UBLAS_INLINE
        reverse_iterator rbegin () {
            return reverse_iterator (end ());
        }
        BOOST_UBLAS_INLINE
        reverse_iterator rend () {
            return reverse_iterator (begin ());
        }

        // Allocator
        allocator_type get_allocator () {
            return alloc_;
        }

    private:
        // Provide destroy as a non member function
        BOOST_UBLAS_INLINE
        static void static_destroy (reference p) {
            (&p) -> ~value_type ();
        }
        ALLOC alloc_;
        size_type capacity_;
        pointer data_;
        size_type size_;
    };


    namespace detail {
#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
        template<class A, class T>
        struct map_traits {
            typedef BOOST_UBLAS_TYPENAME A::mapped_type &reference;
        };
        template<class I, class T, class ALLOC>
        struct map_traits<map_array<I, T, ALLOC>, T > {
            typedef typename map_array<I, T, ALLOC>::data_reference reference;
        };
#else
#if defined (BOOST_UBLAS_STRICT_MAP_ARRAY)
#error BOOST_UBLAS_STRICT_MAP_ARRAY require partial template speciazation
#endif
        // ISSUE: T is actually only required for VC6 as it can't find mapped_type
        template<class A, class T>
        struct map_traits {
            typedef T &reference;
        };
#endif

        // reserve helpers for map_array and generic maps
        // ISSUE should be in map_traits but want to use on all compilers

        template<class M>
        BOOST_UBLAS_INLINE
        void map_reserve (M &/* m */, typename M::size_type /* capacity */) {
        }
        template<class I, class T, class ALLOC>
        BOOST_UBLAS_INLINE
        void map_reserve (map_array<I, T, ALLOC> &m, typename map_array<I, T, ALLOC>::size_type capacity) {
            m.reserve (capacity);
        }

        template<class M>
        BOOST_UBLAS_INLINE
        typename M::size_type map_capacity (M &m) {
            return m.size ();
        }
        template<class I, class T, class ALLOC>
        BOOST_UBLAS_INLINE
        typename map_array<I, T, ALLOC>::size_type map_capacity (map_array<I, T, ALLOC> &m) {
            return m.capacity ();
        }
    }

    // This specialization is missing in Dinkumware's STL?!
    template<class I, class T, class F, class ALLOC>
    BOOST_UBLAS_INLINE
    void swap (std::map<I, T, F, ALLOC> &a1, std::map<I, T, F, ALLOC> &a2) {
        if (&a1 != &a2)
            a1.swap (a2);
    }


#ifdef BOOST_UBLAS_DEPRACATED
// Depracated due to:
//  no allocator implementation
//  inconsitent value_type zero init
//  non STL typedefs

    // Set array
    template<class I, class ALLOC>
    class set_array {
    public:
        typedef ALLOC allocator_type;
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef I index_type;
        typedef I value_type;
        typedef const I &const_reference;
        typedef I &reference;
        typedef const I *const_pointer;
        typedef I *pointer;

        // Construction and destruction
        BOOST_UBLAS_INLINE
        set_array ():
            capacity_ (0), data_ (new value_type [0]), size_ (0) {
        }
        BOOST_UBLAS_INLINE
        set_array (const set_array &a):
            capacity_ (a.size_), data_ (new value_type [a.size_]), size_ (a.size_) {
            *this = a;
        }
        BOOST_UBLAS_INLINE
        ~set_array () {
            delete [] data_;
        }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size) {
            BOOST_UBLAS_CHECK (size_ <= capacity_, internal_logic ());
            if (size > capacity_) {
                pointer data = new value_type [size << 1];
                std::copy (data_, data_ + (std::min) (size, size_), data);
                std::fill (data + (std::min) (size, size_), data + size, value_type (0));
                delete [] data_;
                capacity_ = size << 1;
                data_ = data;
            }
            size_ = size;
            BOOST_UBLAS_CHECK (size_ <= capacity_, internal_logic ());
        }

        // Reserving
        BOOST_UBLAS_INLINE
        void reserve (size_type capacity) {
            BOOST_UBLAS_CHECK (size_ <= capacity_, internal_logic ());
            if (capacity > capacity_) {
                pointer data = new value_type [capacity];
                std::copy (data_, data_ + size_, data);
                delete [] data_;
                capacity_ = capacity;
                data_ = data;
            }
            BOOST_UBLAS_CHECK (size_ <= capacity_, internal_logic ());
        }

        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator [] (index_type i) {
            pointer it = find (i);
            if (it == end ())
                it = insert (end (), i);
            BOOST_UBLAS_CHECK (it != end (), internal_logic ());
            return *it;
        }

        // Assignment
        BOOST_UBLAS_INLINE
        set_array &operator = (const set_array &a) {
            if (this != &a) {
                resize (a.size_);
                std::copy (a.data_, a.data_ + a.size_, data_);
            }
            return *this;
        }
        BOOST_UBLAS_INLINE
        set_array &assign_temporary (set_array &a) {
            swap (a);
            return *this;
        }

        // Swapping
        BOOST_UBLAS_INLINE
        void swap (set_array &a) {
            if (this != &a) {
                std::swap (capacity_, a.capacity_);
                std::swap (data_, a.data_);
                std::swap (size_, a.size_);
            }
        }
#ifndef BOOST_UBLAS_NO_MEMBER_FRIENDS
        BOOST_UBLAS_INLINE
        friend void swap (set_array &a1, set_array &a2) {
            a1.swap (a2);
        }
#endif

        // Element insertion and deletion
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer push_back (pointer it, const value_type &p) {
            if (size () == 0 || (*(it = end () - 1)) < p) {
                resize (size () + 1);
                *(it = end () - 1) = p;
                return it;
            }
            external_logic ().raise ();
            return it;
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer insert (pointer it, const value_type &p) {
            it = detail::lower_bound (begin (), end (), p, std::less<value_type> ());
            difference_type n = it - begin ();
            BOOST_UBLAS_CHECK (size () == 0 || size () == size_type (n) || *it != p, external_logic ());
            resize (size () + 1);
            it = begin () + n;
            std::copy_backward (it, end () - 1, end ());
            *it = p;
            return it;
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        void erase (pointer it) {
            BOOST_UBLAS_CHECK (begin () <= it && it < end (), bad_index ());
            std::copy (it + 1, end (), it);
            resize (size () - 1);
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        void erase (pointer it1, pointer it2) {
            BOOST_UBLAS_CHECK (begin () <= it1 && it1 < it2 && it2 <= end (), bad_index ());
            std::copy (it2, end (), it1);
            resize (size () - (it2 - it1));
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        void clear () {
            resize (0);
        }

        // Element lookup
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_pointer find (index_type i) const {
            const_pointer it (detail::lower_bound (begin (), end (), i, std::less<value_type> ()));
            if (it == end () || *it != i)
                it = end ();
            return it;
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer find (index_type i) {
            pointer it (detail::lower_bound (begin (), end (), i, std::less<value_type> ()));
            if (it == end () || *it != i)
                it = end ();
            return it;
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        const_pointer lower_bound (index_type i) const {
            return detail::lower_bound (begin (), end (), i, std::less<value_type> ());
        }
        // This function seems to be big. So we do not let the compiler inline it.
        // BOOST_UBLAS_INLINE
        pointer lower_bound (index_type i) {
            return detail::lower_bound (begin (), end (), i, std::less<value_type> ());
        }

        // Iterators simply are pointers.

        typedef const_pointer const_iterator;

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return data_;
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return data_ + size_;
        }

        typedef pointer iterator;

        BOOST_UBLAS_INLINE
        iterator begin () {
            return data_;
        }
        BOOST_UBLAS_INLINE
        iterator end () {
            return data_ + size_;
        }

        // Reverse iterators

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef std::reverse_iterator<const_iterator, value_type, const_reference> const_reverse_iterator;
#else
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
#endif

        BOOST_UBLAS_INLINE
        const_reverse_iterator rbegin () const {
            return const_reverse_iterator (end ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator rend () const {
            return const_reverse_iterator (begin ());
        }

#ifdef BOOST_MSVC_STD_ITERATOR
        typedef std::reverse_iterator<iterator, value_type, reference> reverse_iterator;
#else
        typedef std::reverse_iterator<iterator> reverse_iterator;
#endif

        BOOST_UBLAS_INLINE
        reverse_iterator rbegin () {
            return reverse_iterator (end ());
        }
        BOOST_UBLAS_INLINE
        reverse_iterator rend () {
            return reverse_iterator (begin ());
        }

        // Allocator
        allocator_type get_allocator () {
            return alloc_;
        }

    private:
        size_type capacity_;
        pointer data_;
        size_type size_;
    };

    // This specialization is missing in Dinkumware's STL?!
    template<class I, class F, class ALLOC>
    BOOST_UBLAS_INLINE
    void swap (std::set<I, F, ALLOC> &a1, std::set<I, F> &a2) {
        if (&a1 != &a2)
            a1.swap (a2);
    }
#endif

}}}

#endif
