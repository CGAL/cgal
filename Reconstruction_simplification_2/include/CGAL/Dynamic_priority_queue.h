#ifndef DYNAMICS_PRIORITY_QUEUE_H_
#define DYNAMICS_PRIORITY_QUEUE_H_

/*
 * class DynamicPriorityQueue
 * USAGE:
 *    + Implement abstract method 'compare(const T&, const T&)'
 *    + Class T must have operators '=' and '<'
 */

#include <map>
#include <vector>
#include <cassert>

template <class T>
class Dynamic_priority_queue {
protected:  
    typedef std::map<T, int> Map;
    std::vector<T> _queue;
    Map _elements;
    
public:
    Dynamic_priority_queue() {

    }
    
    virtual ~Dynamic_priority_queue() {
        clean();
    }
    
    bool empty() const {
        return _queue.empty();
    }
    
    unsigned int size() const {
        return _queue.size();
    }
    
    const T& get_element(int i) const {
        assert(i >= 0);
        assert(i < size());
        return _queue[i];
    }
    
    bool contains(const T& x) {
        return (get_position(x) >= 0);
    }
    
    void clean() {
        _queue.clear();
        _elements.clear();
    }
    
    void push(const T& x) {
        int pos = insert(x);
        move_up(pos);
    }
    
    bool update_if_better(const T& x) {
        int pos = get_position(x);
        if (pos < 0) return false;
        
        const T& old = get_element(pos);
        if (compare(old, x)) return false;
        
        return update(x);
    }
    
    bool update(const T& x) {
        bool ok = remove(x);
        if (ok) push(x);
        return ok;
    }
    
    bool remove(const T& x) {
        int pos = get_position(x);
        if (pos < 0) return false;
        
        T old = get_element(pos);
        erase(pos);
        
        if (size() < 2) return true;
        if (pos == size()) return true;
        
        if (compare(old, get_element(pos)))
            move_down(pos);
        else 
            move_up(pos);
        return true;
    }
    
    bool extract(T& x) {
        if (empty()) return false;
        x = top();
        pop();
        return true;
    }
    
    T top() const {
        return _queue.front();
    }
    
    void pop() {
        if (empty()) return;
        erase(0);
        if (empty()) return;
        move_down(0);
    }
    
    virtual bool compare(const T& a, const T& b) const = 0;
    
protected:
    int get_position(const T& x) {
        typename Map::const_iterator it = _elements.find(x);
        if (it == _elements.end()) return -1;
        return it->second;
    }
    
    int insert(const T& x) {
        int pos = int( size() );
        _queue.push_back(x);
        _elements[x] = pos;
        return pos;
    }
    
    void erase(int pos) {
        swap(pos, size()-1);
        drop_last_element();
    }
    
    void drop_last_element() {
        T last = _queue.back();
        _queue.pop_back();
        _elements.erase(last);
    }
    
    void swap(int i, int j) {
        T x_i = get_element(i);
        T x_j = get_element(j);
        place(x_i, j);
        place(x_j, i);
    }
    
    void place(const T& x, int i) {
        _queue[i] = x;
        _elements[x] = i;
    }
    
    int parent(int i) const {
        if (i%2 == 0) return (i-2)/2;
        return (i-1)/2;
    }
    
    int left(int i) const {
        return 2*i + 1;
    }
    
    int right(int i) const {
        return 2*i + 2;
    }
    
    void move_down(int i) {
        int pos = i;
        int l = left(pos);
        int r = right(pos);
        T moving = get_element(pos);
        
        int better;
        int N = int(size());
        while (l < N) {
            if ( r < N && compare(get_element(r), get_element(l)) )
                better = r;
            else
                better = l;
            
            if ( compare(moving, get_element(better)) ) break;
            
            place(get_element(better), pos);
            pos = better;
            l = left(pos);
            r = right(pos);
        }
        if (pos != i) place(moving, pos);
    }
    
    void move_up(int i) {
        int pos = i;
        int p = parent(pos);
        T moving = get_element(pos);
        
        while (pos > 0 && compare(moving, get_element(p))) {
            place(get_element(p), pos);
            pos = p;
            p = parent(pos);
        }
        if (pos != i) place(moving, pos);
    }

};


#endif
