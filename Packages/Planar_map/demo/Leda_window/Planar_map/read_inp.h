#ifndef READ_INP_H
#define READ_INP_H

#include <CGAL/basic.h>
#include <CGAL/assertions.h>

template <class Traits>
class PM_input
{
public:
  typedef typename Traits::Point_2 Point;
  typedef struct { int s, t; } Curve;

  PM_input()
  {
    num_pnts = num_cvs = 0;
    pnts = NULL;
    cvs = NULL;
  }

  ~PM_input()
  {
    delete_all();
  }

  void delete_all()
  {
    if (pnts != NULL)
      delete[] pnts;
    if (cvs != NULL)
      delete[] cvs;
    num_pnts = num_cvs = 0;
    pnts = NULL;
    cvs = NULL;
  }

  void alloc(int np, int nc)
  {
    delete_all();
    pnts = new Point[np];
    num_pnts = np;
    cvs = new Curve[nc];
    num_cvs = nc;
  }

  void set(int i, int source, int target)
  {
    if ((i < 0) || (i >= num_cvs))
      return;
    cvs[i].s = source;
    cvs[i].t = target;
  }

  void set(int i, const Point &p)
  {
    if ((i < 0) || (i >= num_pnts))
      return;
    pnts[i] = p;
  }

  bool get(int i, int &source, int &target) const
  {
    if ((i < 0) || (i >= num_cvs))
      return false;
    source = cvs[i].s;
    target = cvs[i].t;
    return true;
  }

  bool get(int i, Point &p) const
  {
    if ((i < 0) || (i >= num_pnts))
      return false;
    p = pnts[i];
    return true;
  }
  
  int curve_source(int i) const
  {
    if ((i < 0) || (i >= num_cvs))
      return -1;
    return cvs[i].s;
  }

  int curve_target(int i) const
  {
    if ((i < 0) || (i >= num_cvs))
      return -1;
    return cvs[i].t;
  }

  int get_num_pnts() const { return num_pnts; }
  
  int get_num_cvs() const { return num_cvs; }
  
  Point point(int i) const
  {
    if ((i < 0) || (i >= num_pnts))
      return Point(0,0);
    return pnts[i];
  }
  
  friend std::ostream &operator<<(std::ostream &os, const PM_input &pmi)
  {
    CGAL::set_ascii_mode(os);
    os << pmi.num_pnts << std::endl;
    int i;
    for (i = 0; i < pmi.num_pnts; i++)
    {
      os << pmi.pnts[i] << std::endl;
    }
    os << pmi.num_cvs << std::endl;
    for (i = 0; i < pmi.num_cvs; i++)
    {
      os << pmi.cvs[i].s << "  " << pmi.cvs[i].t << std::endl;
    }
    return os;
  }

  friend std::istream &operator>>(std::istream &is, PM_input &pmi)
  {
    CGAL::set_ascii_mode(is);
    pmi.delete_all();
    Point p;
    
    is >> pmi.num_pnts;
#ifdef CGAL_PM_READ_DEBUG
    std::cerr << "pmi.num_pnts " << pmi.num_pnts << std::endl;
#endif
    pmi.pnts = new Point[pmi.num_pnts];
    CGAL_assertion(pmi.pnts != NULL);
    int i;
#ifdef USE_LEDA_RAT_KERNEL
    number_type x,y;
#endif
    for (i = 0; i < pmi.num_pnts; i++)
    {
#ifndef USE_LEDA_RAT_KERNEL
      is >> pmi.pnts[i];
#else
      is >> x >> y;
      pmi.pnts[i]=Point(x,y);
#endif
    }
		
    is >> pmi.num_cvs;
#ifdef CGAL_PM_READ_DEBUG
    std::cerr << "pmi.num_cvs " << pmi.num_cvs << std::endl;
#endif
    pmi.cvs = new Curve[pmi.num_cvs];
    CGAL_assertion(pmi.cvs != NULL);
    for (i = 0; i < pmi.num_cvs; i++)
    {
      is >> pmi.cvs[i].s >> pmi.cvs[i].t;
    }
    
    return is;
  }

private:
  Point *pnts;
  Curve *cvs;
  int num_pnts;
  int num_cvs;
};

#endif
