/*  in_list_prog.C                  */
/*  ------------------------------- */
#include <CGAL/basic.h>
#include <assert.h>
#include <algobase.h>
#include <CGAL/In_place_list.h>

struct item : public CGAL_In_place_list_base<item> {
    int key;
    item() {}
    item( const item& i) : CGAL_In_place_list_base<item>(i), key(i.key) {}
    item( int i) : key(i) {}
    bool operator== (const item& i) const { return key == i.key;}
    bool operator!= (const item& i) const { return ! (*this == i);}
    bool operator== (int i) const         { return key == i;}
    bool operator!= (int i) const         { return ! (*this == i);}
    bool operator<  (const item& i) const { return key < i.key;}
};

main() {
    typedef CGAL_In_place_list<item,true> List;
    List  l;
    item* p = new item(1);
    l.push_back(  *p);
    l.push_back(  *new item(2));
    l.push_front( *new item(3));
    l.push_front( *new item(4));
    l.push_front( *new item(2));
    List::iterator i = l.begin();
    ++i;
    l.insert(i,*new item(5));
    l.insert(p,*new item(5));
    int a[7] = {2,5,4,3,5,1,2};
    assert( equal( l.begin(), l.end(), a));
    #ifndef CGAL_CFG_NO_LAZY_INSTANTIATION
    l.sort();
    l.unique();
    int b[5] = {1,2,3,4,5};
    assert( l.size() == 5);
    assert( equal( l.begin(), l.end(), b));
    #endif
    return 0;
}
