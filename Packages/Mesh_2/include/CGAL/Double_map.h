#include <map>
#include <pair>

template <class _Key, class _Data, class _Direct_order = std::less<Key>, 
	  class _Reverse_order = std::less<Data> >
class Double_map
{
public:
  typedef _Key Key;
  typedef _Data Data;
  typedef _Direct_order Direct_order;
  typedef _Reverse_order Reverse_order;
  
  typedef std::map <Key, Data, Direct_order> Direct_map;
  typedef std::multimap <Data, Key, Reverse_order> Reverse_multi_map;

  typedef std::pair<Key, Data> Direct_entry;
  typedef std::pair<Data, Key> Reverse_entry;

  typedef typename Direct_entry::size_type size_type;

  typedef typename Reverse_multi_map::iterator reverse_iterator;
  typedef typename Direct_map::iterator direct_iterator;
  typedef reverse_iterator iterator;

private:
  // Private members
  Direct_map direct_func;
  Reverse_multi_map reverse_func;

public :
  // The default constructor
  Double_map () {}


  // Queries
  bool is_empty() const
  {
      return(direct_func.empty());
  }

  size_type size() const
  {
    CGAL_assertion(is_valid());
    return direct_func.size();
  }

  bool is_valid() const
  {
    return(direct_func.size()==reverse_func.size());
  }
  
  void clear()
    {
      direct_func.clear();
      reverse_func.clear();
    }

  // Assignation
  bool insert(Key& k, Data d)
    {
      std::pair<direct_iterator, bool> result;

      result=direct_func.insert(Direct_entry (k, d));

      if (result.second==true)
	reverse_func.insert(Reverse_entry (d,k));

      return(result.second);
    }

  bool erase(Key& k)
    {
      std::pair<reverse_iterator, reverse_iterator> eq_keys;
      reverse_iterator mit;

      direct_iterator pos=direct_func.find(k);
      if (pos!=direct_func.end())
	{
	  Data d=direct_func[k];
	  direct_func.erase(pos);

	  eq_keys=reverse_func.equal_range(d);
	  mit=eq_keys.first;
	  while(mit!=eq_keys.second)
	    {
	      reverse_iterator temp_mit;
	      mit++;
	      temp_mit=mit;
	      mit--;

	      if (mit->second==k)
		reverse_func.erase(mit);

	      mit=temp_mit;
	    }
	  if (mit->second==k)
	    reverse_func.erase(mit);
	}

      return(false);
    }

  // Access
  reverse_iterator front() const
    {
      return(reverse_func.begin());
    }
};
