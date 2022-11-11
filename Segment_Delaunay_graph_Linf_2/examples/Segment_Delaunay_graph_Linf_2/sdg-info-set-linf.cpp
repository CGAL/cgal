// example that shows how to add info to input sites and how this is
// propagated using the storage traits with info
//
// the input sites are considered to have a label (a character
// associated with them)
// points on the plane belonging to sites of different labels get all
// the labels from the different sites they belong to.

// standard includes
#include <iostream>
#include <string>
#include <fstream>
#include <list>
#include <cassert>

// a class representing a set of info items
template<class Info_item>  struct Info_set_merge_info;

template<class Info_item_t>
class Info_set
{
public:
  typedef Info_item_t                         Info_item;

private:
  friend struct Info_set_merge_info<Info_item>;

  typedef std::list<Info_item>                Info_list;

public:
  typedef typename Info_list::const_iterator  Info_item_iterator;
  typedef typename Info_list::size_type       size_type;

  Info_set() {}

  Info_set(Info_item info) {
    info_list_.push_back(info);
  }

  template<typename InputIterator>
  Info_set(InputIterator first, InputIterator beyond)
    : info_list_(first, beyond) {}

  size_type size() const { return info_list_.size(); }
  bool is_empty() const { return info_list_.empty(); }

  Info_item_iterator info_items_begin() const {
    return info_list_.begin();
  }

  Info_item_iterator info_items_end() const {
    return info_list_.end();
  }

private:
  // private constructor from list of info items
  Info_set(Info_list info_list) : info_list_(info_list) {}

  // private access to list of info items
  const Info_list& info_list() const { return info_list_; }

private:
  Info_list info_list_;
};


// output operator for the set of info items; it assumes that the
// output operator is defined for info items
template<class Info_item>
std::ostream&
operator<<(std::ostream& os, const Info_set<Info_item>& info)
{
  if ( info.is_empty() ) {
    return os << "{}";
  }

  typedef typename Info_set<Info_item>::Info_item_iterator iterator;
  iterator last = --info.info_items_end();
  os << "{";
  for (iterator it = info.info_items_begin(); it != last; ++it) {
    os << *it << ", ";
  }
  os << *last << "}";

  return os;
}

// functor that defines how to convert color info when:
// 1. constructing the storage site of an endpoint of a segment
// 2. a segment site is split into two sub-segments
template<class Info_item_t>
struct Info_set_convert_info
{
  typedef Info_item_t                  Info_item;
  typedef const Info_set<Info_item>&   result_type;

  inline
  result_type operator()(const Info_set<Info_item>& info0, bool) const
  {
    // just return the info of the supporting segment
    return info0;
  }

  inline
  result_type operator()(const Info_set<Info_item>& info0,
                         const Info_set<Info_item>& , bool) const
  {
    // just return the info of the supporting segment
    return info0;
  }
};


// functor that defines how to merge info items when a site (either
// point or segment) corresponds to point(s) on plane belonging to
// more than one input site
template<class Info_item_t>
struct Info_set_merge_info
{
  typedef Info_item_t          Info_item;
  typedef Info_set<Info_item>  result_type;

  inline
  Info_set<Info_item> operator()(const Info_set<Info_item>& info0,
                                 const Info_set<Info_item>& info1) const
  {
    typedef typename Info_set<Info_item>::Info_list       Info_list;

    // return as new info the union of the two infos
    Info_list info_union = info0.info_list();
    Info_list copy = info1.info_list();

    info_union.splice(info_union.end(), copy);
    return info_union;
  }
};


// finally a class that generates info when the info items are
// std::strings
struct Generate_info
{
  static unsigned int ctr;

  template<class Site>
  std::string operator()(const Site& t) const
  {
    if ( t.is_point() ) {
      char c = 'A' + ctr++;
      return std::string(1, c);
    }
    char c1 = 'A' + ctr++;
    char c2 = 'A' + ctr++;
    return std::string(1, c1) + std::string(1, c2);
  }
};

unsigned int Generate_info::ctr = 0;


// choose the kernel
#include <CGAL/Simple_cartesian.h>

struct Rep : public CGAL::Simple_cartesian<double> {};

// typedefs for the geometric traits, storage traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_Linf_hierarchy_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_traits_with_info_2.h>


typedef CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<Rep> Gt;

// define the info and the convert and merge functors
typedef std::string                         Info_item;
typedef Info_set<Info_item>                 Info;
typedef Info_set_convert_info<Info_item>    Convert_info;
typedef Info_set_merge_info<Info_item>      Merge_info;


// define the storage traits with info

typedef
CGAL::Segment_Delaunay_graph_storage_traits_with_info_2<Gt,
                                                             Info,
                                                             Convert_info,
                                                             Merge_info>
ST;

typedef CGAL::Segment_Delaunay_graph_Linf_hierarchy_2<Gt,ST> SDG2;

typedef SDG2::Finite_vertices_iterator     FVIT;
typedef SDG2::Site_2                       Site_2;


int main( int argc, char *argv[] )
{
  if ( ! (( argc == 1 ) || (argc == 2)) ) {
    std::cout <<"usage: "<< argv[0] <<" [filename]\n";
  }

  std::ifstream ifs( (argc == 1) ? "data/sitesxx.cin" : argv[1] );
  assert( ifs );

  SDG2 sdg;
  Site_2 site;
  Generate_info generate;

  // read the sites and their info and insert them in the
  // segment Delaunay graph; print them as you read them
  std::cout << std::endl;
  std::cout << "Input sites:" << std::endl;
  std::cout << "------------" << std::endl;
  while ( ifs >> site ) {
    Info info = generate(site);
    std::cout << site << std::flush;
    std::cout << "\r\t\t\t" << info << std::endl;
    sdg.insert(site, info);
  }
  std::cout << std::endl;

  // validate the segment Delaunay graph
  assert( sdg.is_valid(true) );

  // print the sites of the segment Delaunay graph and their info
  std::cout << std::endl;
  std::cout << "Output sites:" << std::endl;
  std::cout << "-------------" << std::endl;
  for (FVIT it = sdg.finite_vertices_begin();
       it != sdg.finite_vertices_end(); ++it) {
    if ( it->site().is_point() ) {
      std::cout << it->site() << std::flush;
      std::cout << "\r\t\t\t" << it->storage_site().info() << std::endl;
    }
  }
  for (FVIT it = sdg.finite_vertices_begin();
       it != sdg.finite_vertices_end(); ++it) {
    if ( it->site().is_segment() ) {
      std::cout << it->site() << std::flush;
      std::cout << "\r\t\t\t" << it->storage_site().info() << std::endl;
    }
  }
  std::cout << std::endl;

  return 0;
}
