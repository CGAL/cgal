// Try Sylvain proposal  for a fully flexible design
// using the rebind trick
// try a general example (not in CGAL), first
#include <list>

struct TDS_Bidon {
 typedef int  Vertex_handle;
};


//template <class TDS = void> 
template <class TDS = TDS_Bidon>
class  Vertex_base 
{
public:
  typedef typename TDS::Vertex_handle            Vertex_handle;

  template < class My_TDS>
  struct TDS_rebind { typedef Vertex_base<My_TDS> Rebound;};
    
  Vertex_base() : _vh(0){}
  Vertex_base(Vertex_handle vh) : _vh(vh) {}
  void set_vh(Vertex_handle vh) { _vh = vh;}
  Vertex_handle get_vh() {return _vh ;}

private:
  Vertex_handle _vh;
};

template<class Vb>
class TDS_2 {
public:
  typedef TDS_2<Vb>  Self;
  typedef typename Vb::template TDS_rebind<Self>::Rebound Vertex;

  typedef std::list<Vertex>             Vertex_container;
  typedef Vertex*                       Vertex_handle;

protected:
  Vertex_container container;  

public:
  Vertex_handle create_vertex(Vertex_handle vh){
    Vertex* vv = new Vertex(vh);
    container.push_back(*vv);
    return vv;
    }
};



//
// A user defined Vertex 
template <class TDS = TDS_Bidon> 
class  My_vertex_base : public Vertex_base<TDS>
{
public:
  typedef Vertex_base<TDS>                     Base;
  typedef typename TDS::Vertex_handle          Vertex_handle;

  template < class My_TDS>
  struct TDS_rebind { typedef My_vertex_base<My_TDS> Rebound;};
    
  My_vertex_base() : Base (){}
  My_vertex_base(Vertex_handle vh) : Base(vh) {}
      
  void set_wahou(Vertex_handle vh) { wahou = vh;}
  Vertex_handle get_wahou() {return wahou ;}

private:
 Vertex_handle wahou;
};


typedef Vertex_base<>                        Vb;
typedef TDS_2<Vb>                            Tds;
typedef My_vertex_base<>                     Myvb;
typedef TDS_2<Myvb>                          Mytds;

typedef Tds::Vertex                          Vertex;
typedef Tds::Vertex_handle                   Vertex_handle;
typedef Mytds::Vertex                        My_vertex;          
typedef Mytds::Vertex_handle                 My_vertex_handle;

int  main()
{
  Tds tds;
  Vertex_handle vh = new Vertex;
  for(int i = 0 ; i < 10 ; i++) {
    vh = tds.create_vertex(vh);
  }

  Mytds mytds;
  My_vertex_handle myvh = new My_vertex;
  for(int i = 0 ; i < 10 ; i++) {
    myvh = mytds.create_vertex(myvh);
    myvh->set_wahou(myvh);
  }
  return 0;
}




