#ifndef CGAL_VECTOR_H
#define CGAL_VECTOR_H

//#include <iostream>

#ifndef __SGI_STL_INTERNAL_VECTOR_H
#  include <stl_vector.h>
#endif

#include <memory>

#define CGAL_NULL_TMPL_ARGS

#ifdef vector
#  undef vector
#endif

// namespace stlport { // replace eventually with 
__STL_BEGIN_NAMESPACE 

template <class T, class _Alloc >
class vector;

template <class T, class _Alloc>
class vector_iterator;

template <class T, class _Alloc>
class vector_const_iterator {
  public:  
    friend class vector<T,_Alloc>;
    friend class vector_iterator<T,_Alloc>;
    typedef typename _Alloc::const_pointer pointer;
    typedef __STLPORT_STD::random_access_iterator_tag iterator_category;
    typedef T value_type;
    typedef typename _Alloc::difference_type difference_type;
    typedef typename _Alloc::const_reference reference;
    typedef vector_const_iterator<T,_Alloc> Self;
  protected:
    pointer the_it;
    explicit vector_const_iterator (pointer it) :the_it(it) {}
  public:
    vector_const_iterator():the_it(0) {}

    reference operator*() const {return *the_it;}
    pointer operator->() const {return the_it;}
    Self & operator++() { ++the_it; return *this;}
    Self operator++(int) { return Self(the_it++);}
    Self & operator--() { --the_it; return *this;}
    Self operator--(int) { return Self(the_it--);}
    Self & operator+=(difference_type n)
		{ the_it+=n; return *this;}
    Self operator+(difference_type n) const
    		{ return Self(the_it+n);}
    Self & operator-=(difference_type n)
		{ the_it-=n; return *this;}
    Self operator-(difference_type n) const
   		{ return Self(the_it-n);}
    difference_type operator-(Self j) const
   		{ return the_it-j.the_it;}
    bool operator<(Self j) const
   		{ return the_it<j.the_it;}
    bool operator>(Self j) const
   		{ return the_it>j.the_it;}
    bool operator<=(Self j) const
   		{ return the_it<=j.the_it;}
    bool operator>=(Self j) const
   		{ return the_it>=j.the_it;}
    reference operator[](difference_type n) const { return the_it[n];}


friend bool operator== CGAL_NULL_TMPL_ARGS
    (const vector_const_iterator<T,_Alloc> &, const vector_const_iterator<T,_Alloc> &);
};

template <class T, class _Alloc>
bool operator==(const __STLPORT_STD::vector_const_iterator<T,_Alloc> &it1,
	const __STLPORT_STD::vector_const_iterator<T,_Alloc> &it2)
{
    return it1.the_it == it2.the_it;
}


template <class T, class _Alloc>
inline bool operator!=(const 
		       __STLPORT_STD::vector_const_iterator<T,_Alloc> &it1,
		       const 
		       __STLPORT_STD::vector_const_iterator<T,_Alloc> &it2)
{
    return !(it1== it2);
}

template <class T, class _Alloc>
class vector_iterator {
  public:
    friend class vector<T,_Alloc>;
    typedef typename _Alloc::pointer pointer;
    typedef __STLPORT_STD::random_access_iterator_tag iterator_category;
    typedef T value_type;
    typedef typename _Alloc::difference_type difference_type;
    typedef typename _Alloc::reference reference;
    typedef vector_iterator<T,_Alloc> Self;
  protected:
    pointer the_it;
    explicit vector_iterator (pointer it) :the_it(it) {}
  public:
    vector_iterator():the_it(0) {}

    operator vector_const_iterator<T,_Alloc>() const
		{ return vector_const_iterator<T,_Alloc>(the_it); }
    reference operator*() const {return *the_it;}
    pointer operator->() const {return the_it;}
    Self & operator++() { ++the_it; return *this;}
    Self operator++(int) { return Self(the_it++);}
    Self & operator--() { --the_it; return *this;}
    Self operator--(int) { return Self(the_it--);}
    Self & operator+=(difference_type n)
		{ the_it+=n; return *this;}
    Self operator+(difference_type n) const { return Self(the_it+n);}
    Self & operator-=(difference_type n)
		{ the_it-=n; return *this;}
    Self operator-(difference_type n) const
   		{ return Self(the_it-n);}
    difference_type operator-(Self j) const
   		{ return the_it-j.the_it;}
    bool operator<(Self j) const
   		{ return the_it<j.the_it;}
    bool operator>(Self j) const
   		{ return the_it>j.the_it;}
    bool operator<=(Self j) const
   		{ return the_it<=j.the_it;}
    bool operator>=(Self j) const
   		{ return the_it>=j.the_it;}
    reference operator[](difference_type n) const { return the_it[n];}

friend bool operator== CGAL_NULL_TMPL_ARGS
    (const vector_iterator<T,_Alloc> &, const vector_iterator<T,_Alloc> &);

};


template <class T, class _Alloc>
bool operator==(const __STLPORT_STD::vector_iterator<T,_Alloc> &it1,
	const __STLPORT_STD::vector_iterator<T,_Alloc> &it2)
{
    return it1.the_it == it2.the_it;
}


template <class T, class _Alloc>
inline bool operator!=(const __STLPORT_STD::vector_iterator<T,_Alloc> &it1,
	const __STLPORT_STD::vector_iterator<T,_Alloc> &it2)
{
    return !(it1== it2);
}

template <class T, class _Alloc = __STLPORT_STD::allocator<T> >
class vector {
  public:
    typedef T value_type;
    typedef typename _Alloc::pointer pointer;
    typedef typename _Alloc::const_pointer const_pointer;
    typedef typename _Alloc::reference reference;
    typedef typename _Alloc::const_reference const_reference;
    typedef typename _Alloc::size_type size_type;
    typedef typename _Alloc::difference_type difference_type;
    typedef vector_const_iterator<T,_Alloc> const_iterator;
    typedef vector_iterator<T,_Alloc> iterator;
    typedef __STLPORT_STD::reverse_iterator<iterator> reverse_iterator;
    typedef __STLPORT_STD::reverse_iterator<const_iterator>
            const_reverse_iterator;
  protected:
    pointer data_;
    size_type size_;
    size_type capacity_;
    _Alloc alloc_;
  public:    
    iterator begin() {return iterator(data_);}
    const_iterator begin() const {return const_iterator(data_);}
    iterator end() {return iterator(data_+size_);}
    const_iterator end() const {return const_iterator(data_+size_);}
    reverse_iterator rbegin() {return reverse_iterator(end());}
    reverse_iterator rend() {return reverse_iterator(begin());}
    const_reverse_iterator rbegin() const {return const_reverse_iterator(end());}
    const_reverse_iterator rend() const {return const_reverse_iterator(begin());}
    
    size_type size() const  {return size_;}
    size_type max_size() const
    	{ return size_type(-1) / sizeof(T); }
    size_type capacity() const {return capacity_;}
    bool empty() const {return size_ == 0;}
    
    reference operator[] (size_type n) { return data_[n];}
    const_reference operator[] (size_type n) const { return data_[n];}
    
    explicit vector(const _Alloc & alloc = _Alloc())
		: alloc_(alloc), size_(0),data_(0),capacity_(0){}
    explicit vector(size_type n, const T& x = T() , const _Alloc &alloc = _Alloc() );
    vector(const vector & v) ;
    template <class InputIterator>
    vector(InputIterator f, InputIterator l, const _Alloc &alloc = _Alloc() )
        : alloc_(alloc), size_(0),data_(0),capacity_(0)
        {
            typedef typename iterator_traits<InputIterator>::iterator_category category;
            construct(f,l,category());
        }
    ~vector();
    vector & operator=(const vector &o) {assign(o.begin(), o.end()); return *this;}

    void swap(vector &o) ;
    void reserve(size_type n);
    reference front() {return *data_;}
    const_reference front() const {return *data_;}
    reference back() {return data_[size_-1];}
    const_reference back() const {return data_[size_-1];}
    void push_back(const T& t);
    void pop_back();
    iterator insert(iterator pos, const T& x);
    template <class InputIterator>
    void insert(iterator pos, InputIterator f, InputIterator l)
    {
        typedef iterator_traits<InputIterator>::iterator_category category;
        insert_it(pos,f,l,category());
    }
    void insert(iterator pos, int n, const T& x) { insert_n(pos,(size_type)n,x);}
    void insert(iterator pos, long n, const T& x) { insert_n(pos,(size_type)n,x);}
    void insert(iterator pos, unsigned int n, const T& x) { insert_n(pos,(size_type)n,x);}
    void insert(iterator pos, unsigned long n, const T& x) { insert_n(pos,(size_type)n,x);}
    iterator erase(iterator pos);
    iterator erase(iterator first, iterator last);
    void clear();
    void resize(size_type n,const T& t = T());
    template <class InputIterator>
    void assign(InputIterator f, InputIterator l)
    {
        typedef iterator_traits<InputIterator>::iterator_category category;
        clear();
        construct(f,l,category());
    }
    void assign(size_type n, const T& x);
  protected:
    void move_negative(pointer source, pointer source_end, pointer destination);
    void move_positive(difference_type shift, pointer &source);
    void insert_n(iterator pos, size_type n, const T& x);
    template <class InputIterator>
    void construct(InputIterator cur, InputIterator beyond, __STLPORT_STD::input_iterator_tag)
    {  for ( ; cur != beyond; ++cur) push_back(*cur); }
    template <class ForwardIterator>
    void construct(ForwardIterator cur, ForwardIterator beyond, __STLPORT_STD::forward_iterator_tag)
    {
    	size_type n = distance(cur,beyond);
    	reserve(n);
    	pointer fill = data_ + size_;
    	for (size_type i = 0; i < n; ++i,++fill, ++cur) {
    	    alloc_.construct(fill,*cur);
    	}
    	size_ = n;
    }
    template <class InputIterator>
    void insert_it(iterator pos, InputIterator f,
                   InputIterator l, __STLPORT_STD::input_iterator_tag)
    {
        for (;f != l; ++f) {
            pos = insert(pos, *f);
            ++pos;
        }   
    }
    template <class InputIterator>
    void insert_it(iterator pos, InputIterator f,
                   InputIterator l, __STLPORT_STD::forward_iterator_tag)
    {
        pointer result = pos.the_it;
        size_type n = distance(f,l);
        move_positive(n, result);
        for (size_type i =0; i < n; ++i, ++result, ++f) {
            alloc_.construct(result, *f);
        }
    }

friend bool operator== CGAL_NULL_TMPL_ARGS
    (const vector<T> &, const vector<T> &);

friend bool operator< CGAL_NULL_TMPL_ARGS
    (const vector&v1, const vector &v2);
};


template <class T>
bool operator==(const __STLPORT_STD::vector<T> &v1,
	const __STLPORT_STD::vector<T> &v2)
{
    if (v1.size() != v2.size())
        return false;
    typedef typename __STLPORT_STD::vector<T>::const_iterator ci;
    ci cur1 = v1.begin(), cur2 = v2.begin(), end1 = v1.end();
    for ( ; cur1 != end1; ++cur1, ++cur2) {
        if ( !(*cur1 == *cur2)) return false;
    }
    return true;
}

template <class T>
bool operator< (const __STLPORT_STD::vector<T>&v1, 
		const __STLPORT_STD::vector<T>& v2)
{
    return __STLPORT_STD::lexicographical_compare(v1.begin(), v1.end(),
						  v2.begin(),v2.end());
}

# ifdef __STL_USE_SEPARATE_RELOPS_NAMESPACE

template <class T>
inline bool 
operator!=(const __STLPORT_STD::vector<T>& __x, 
	   const __STLPORT_STD::vector<T>& __y)
{
  return !(__x == __y);
}

template <class T>
inline bool operator>(const __STLPORT_STD::vector<T>& __x, 
		      const __STLPORT_STD::vector<T>& __y)
{
  return __y < __x;
}

template <class T>
inline bool operator<=(const __STLPORT_STD::vector<T>& __x, 
		       const __STLPORT_STD::vector<T>& __y)
{
  return !(__y < __x);
}

template <class T>
inline bool operator>=(const __STLPORT_STD::vector<T>& __x, 
		       const __STLPORT_STD::vector<T>& __y)
{
  return !(__x < __y);
}

#endif

// } // namespace stlport ends ; replace eventually with 
  // 

template <class T, class _Alloc>
vector<T,_Alloc>::
vector(size_type n, const T& x , const _Alloc &alloc)
    : alloc_(alloc)
{
    if ( n == 0) {
	size_ = capacity_ = 0;
	data_ = 0;
    } else {
    	size_ = capacity_ = n;
    	data_ = alloc_.allocate(capacity_);
    	T *p = data_;
    	for (size_type i=0; i<capacity_; ++i, ++p) {
	    alloc_.construct(p,x);
    	}
    }
}

template <class T, class _Alloc>
vector<T,_Alloc>::
vector(const vector & v)
: alloc_(v.alloc_)
{
    if (v.size_ == 0) {
	size_ = capacity_ = 0;
	data_ = 0;
    } else {
	size_ = capacity_ = v.size_;
    	data_ = alloc_.allocate(capacity_);
        pointer p = data_;
        for (size_type i=0; i<capacity_; ++i, ++p) {
	    alloc_.construct(p,v[i]);
        }
    }
}

template <class T, class _Alloc>
vector<T,_Alloc>::
    ~vector()
{
    if (data_ != 0) {
	T *p = data_;
        for (size_type i=0; i<size_; ++i, ++p) {
	    alloc_.destroy(p);
        }
	alloc_.deallocate(data_, capacity_);
    }
}

template <class T, class _Alloc>
void vector<T,_Alloc>::
swap(vector &o)
{
   pointer tmp_data= o.data_;
   o.data_ = data_;
   data_ = tmp_data;
   size_t tmp_size = o.size_;
   o.size_ = size_;
   size_ = tmp_size;
   size_t tmp_capacity = o.capacity_;
   o.capacity_ = capacity_;
   capacity_ = tmp_capacity;
}

template <class T, class _Alloc>
void vector<T,_Alloc>::
reserve(size_type n)
{
    if (n <= capacity_) return;
    pointer new_data = alloc_.allocate(n);
    if (data_ != 0) {
    	pointer old_data = data_;
    	pointer new_data_p = new_data;
    	while (old_data < data_+size_) {
	    alloc_.construct(new_data_p, *old_data);
	    alloc_.destroy(old_data);
	    ++old_data; ++new_data_p;
    	}
    	alloc_.deallocate(data_,capacity_);
    }
    data_ = new_data;
    capacity_ = n;
}

template <class T, class _Alloc>
void vector<T,_Alloc>::
push_back(const T& t)
{
    if (size_ == capacity_) {
	size_t new_capacity = 1.3 * capacity_;
	if (new_capacity < capacity_+10)
    	    new_capacity = capacity_+10;	   
	reserve(new_capacity);
    }
    alloc_.construct(data_+size_,t);
    ++size_;
}

template <class T, class _Alloc>
void vector<T,_Alloc>::
pop_back()
{
    --size_;
    alloc_.destroy(data_ + size_);
}


template <class T, class _Alloc>
void vector<T,_Alloc>::
move_positive(difference_type shift, pointer &source)
{
    pointer result = source;
    size_type new_size = size_ + shift;
    if (new_size > capacity_) {
        pointer new_data = alloc_.allocate(new_size);
        if (data_ != 0) {
    	    pointer from = data_;
    	    pointer dest = new_data;
    	    while (from < source) {
	            alloc_.construct(dest, *from);
	            alloc_.destroy(from);
	            ++from; ++dest;
    	    }
    	    result = dest;
    	    dest += shift;
    	    while (from < data_+size_) {
	            alloc_.construct(dest, *from);
	            alloc_.destroy(from);
	            ++from; ++dest;
    	    }
    	    alloc_.deallocate(data_,capacity_);
        } else {
            result = new_data;
        }
        data_ = new_data;
        capacity_ = new_size;
    } else {
        pointer from = data_ + size_; // 1 past the end
        pointer dest = from + shift;
        while (from != source) {
            --from; --dest;
            alloc_.construct(dest, *from);
            alloc_.destroy(from);
        }
    }
    size_ = new_size;
    source = result;
}

template <class T, class _Alloc>
vector_iterator<T,_Alloc> vector<T,_Alloc>::
insert(iterator pos, const T& x)
{
    pointer result = pos.the_it;
    move_positive(1, result);
    alloc_.construct(result, x);
    return static_cast<iterator>(result);
}

template <class T, class _Alloc>
void vector<T,_Alloc>::
insert_n(iterator pos, size_type n, const T& x)
{
    pointer result = pos.the_it;
    move_positive(n, result);
    for (size_type i =0; i < n; ++i, ++result) {
        alloc_.construct(result, x);
    }
    
}

template <class T, class _Alloc>
void vector<T,_Alloc>::
move_negative(pointer source, pointer source_end, pointer destination)
{
    for ( ; source != source_end; ++source, ++ destination) {
        alloc_.construct(destination, *source);
        alloc_.destroy(source);
    }
}


template <class T, class _Alloc>
vector_iterator<T,_Alloc> vector<T,_Alloc>::
erase(iterator pos)
{
    pointer curpos = pos.the_it;
    alloc_.destroy(curpos);
    move_negative(curpos+1, data_+size_,curpos);
    --size_;
    return pos;
}
    
template <class T, class _Alloc>
vector_iterator<T,_Alloc> vector<T,_Alloc>::
erase(iterator first, iterator last)
{
    pointer firstpos = first.the_it;
    pointer lastpos = last.the_it;
    difference_type diff = lastpos-firstpos;
    for (pointer curpos=firstpos; curpos != lastpos; ++curpos)
        alloc_.destroy(curpos);
    move_negative(firstpos+diff, data_+size_,firstpos);
    size_ -= diff;
    return firstpos;
}


template <class T, class _Alloc>
void vector<T,_Alloc>::
clear()
{
    if (size_ == 0) return;
    pointer cur = data_;
    for (size_type i = 0; i < size_; ++i, ++cur) {
        alloc_.destroy(cur);
    }
    size_ = 0;
}

template <class T, class _Alloc>
void vector<T,_Alloc>::
resize(size_type n,const T& t)
{
    if (n == size_) return;
    if (n < size_) {
        pointer cur = data_ + n;
        for ( ; n < size_; --size_, ++cur) {
            alloc_.destroy(cur);
        }
    } else {
        reserve(n);
        while (size_ < n)
            push_back(t);
    }
}

template <class T, class _Alloc>
void vector<T,_Alloc>::
assign(size_type n, const T& x)
{
    clear();
    reserve(n);
    for (size_t i = 0; i < n; ++i)
        push_back(x);
}
__STL_END_NAMESPACE

#endif
