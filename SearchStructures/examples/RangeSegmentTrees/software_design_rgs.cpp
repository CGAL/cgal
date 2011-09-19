#include <CGAL/Segment_tree_d.h>
#include <CGAL/Tree_traits.h>

using namespace CGAL;

struct Data{
  int min,max;
  double point;
};

struct Window{
  int min,max;
  double min_point, max_point;
};

class Interval_traits{
 public:
  typedef int Key;
  Key get_left( const Data&  d) const {return d.min;}
  Key get_right( const Data&  d) const {return d.max;}
  Key get_left_win( const Window& w) const {return w.min;}
  Key get_right_win(  const Window& w) const {return w.max;}
  static bool comp( const Key& key1, const Key& key2) {return (key1 < key2);}
};

typedef Tree_anchor<Data,Window> Tree_Anchor;
typedef Segment_tree_d<Data,Window,Interval_traits> Segment_Tree_d;

void remove_warning(Segment_Tree_d*){}

int main(){
  Tree_Anchor *anchor = new Tree_Anchor;
  Segment_Tree_d *segment_tree = new Segment_Tree_d(*anchor);
  remove_warning(segment_tree);
}
