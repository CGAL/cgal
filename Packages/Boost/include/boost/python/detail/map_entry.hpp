// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef MAP_ENTRY_DWA2002118_HPP
# define MAP_ENTRY_DWA2002118_HPP

namespace boost { namespace python { namespace detail { 

// A trivial type that works well as the value_type of associative
// vector maps
template <class Key, class Value>
struct map_entry
{
    map_entry() {}
    map_entry(Key k) : key(k), value() {}
    map_entry(Key k, Value v) : key(k), value(v) {}
    
    bool operator<(map_entry const& rhs) const
    {
        return this->key < rhs.key;
    }
        
    Key key;
    Value value;
};

template <class Key, class Value>
bool operator<(map_entry<Key,Value> const& e, Key const& k)
{
    return e.key < k;
}

template <class Key, class Value>
bool operator<(Key const& k, map_entry<Key,Value> const& e)
{
    return k < e.key;
}


}}} // namespace boost::python::detail

#endif // MAP_ENTRY_DWA2002118_HPP
