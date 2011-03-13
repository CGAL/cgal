#include <cassert>
#include <algorithm>
#include <CGAL/In_place_list.h>

using CGAL::In_place_list_base;

struct item : public In_place_list_base<item> {
    int key;
    item() {}
    item( const item& i) : In_place_list_base<item>(i), key(i.key) {}
    item( int i) : key(i) {}
    bool operator== (const item& i) const { return key == i.key;}
    bool operator!= (const item& i) const { return ! (*this == i);}
    bool operator== (int i) const         { return key == i;}
    bool operator!= (int i) const         { return ! (*this == i);}
    bool operator<  (const item& i) const { return key < i.key;}
};

int main() {
    typedef CGAL::In_place_list<item,true> List;
    List l;
    item* p = new item(1);
    l.push_back(*p);
    l.push_back(*new item(2));
    l.push_front(*new item(3));
    l.push_front(*new item(4));
    l.push_front(*new item(2));
    List::iterator i = l.begin();
    ++i;
    l.insert(i, *new item(5));
    l.insert(p, *new item(5));
    int a[7] = {2,5,4,3,5,1,2};
    bool ok = std::equal(l.begin(), l.end(), a);
    assert(ok);
    l.sort();
    l.unique();
    assert(l.size() == 5);
    int b[5] = {1,2,3,4,5};
    ok = std::equal(l.begin(), l.end(), b);
    assert(ok);
    return 0;
}
