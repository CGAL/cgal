// Copyright (c) 2007  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@gmail.com>

#if defined(_MSVC_VER)
#  pragma warning(disable:4503) // decorated name length exceeded, name was truncated
#endif

#include "basics.h"
#include "test_self_intersection.h"

#include <memory>

//#define TEST_TEST_TRACE_ENABLED

template<std::size_t N> struct eat_sizeof {};

#ifdef TEST_TEST_TRACE_ENABLED
#  define TEST_TRACE(m) std::cerr << m << std::endl;
#else
#  define TEST_TRACE(m) (eat_sizeof<sizeof(m)>())

#endif


// This is here only to allow a breakpoint to be placed so I can trace back the problem.
void error_handler (char const* what, char const* expr, char const* file, int line, char const* msg)
{
  cerr << "CGAL error: " << what << " violation!" << endl
       << "Expr: " << expr << endl
       << "File: " << file << endl
       << "Line: " << line << endl;
  if(msg != nullptr)
    cerr << "Explanation:" << msg << endl;

  throw std::runtime_error(expr);
}

namespace SMS = CGAL::Surface_mesh_simplification;

typedef SMS::Edge_profile<Surface> Profile;

typedef std::shared_ptr<Surface> SurfaceSP;

// Constructs a flat polyhedron containing just the link of an edge or vertex.
class Link_builder : public CGAL::Modifier_base<Surface::HalfedgeDS>
{
  map<std::size_t, std::size_t> mMap;

  int mVIdx;

protected:

  typedef boost::graph_traits<Surface> GraphTraits;

  typedef CGAL::Halfedge_around_source_iterator<Surface> out_edge_iterator;
  typedef GraphTraits::halfedge_descriptor   halfedge_descriptor;
  typedef GraphTraits::vertex_descriptor vertex_descriptor;
  typedef CGAL::Polyhedron_incremental_builder_3<Surface::HalfedgeDS> Builder;

  Link_builder() : mVIdx(0) {}

  void add_vertex(Builder& aBuilder, Vertex_const_handle v)
  {
    aBuilder.add_vertex(v->point());
    mMap.insert(make_pair(v->id(),mVIdx++));
  }

  void add_triangle(Builder& aBuilder,
                    const Profile::Triangle& aTri)
  {
    aBuilder.begin_facet();
    aBuilder.add_vertex_to_facet(mMap[aTri.v0->id()]);
    aBuilder.add_vertex_to_facet(mMap[aTri.v1->id()]);
    aBuilder.add_vertex_to_facet(mMap[aTri.v2->id()]);
    aBuilder.end_facet();
  }
};

// Constructs a flat polyhedron containing just the link of an edge
class Edge_link_builder : public Link_builder
{
  const Profile& mProfile;


public:

  Edge_link_builder(const Profile& aProfile) : mProfile(aProfile) {}

  void operator()(Surface::HalfedgeDS& hds)
  {
    Builder B(hds, true);

    B.begin_surface(2 + mProfile.link().size(), mProfile.triangles().size());

    this->add_vertex(B,mProfile.v0());
    this->add_vertex(B,mProfile.v1());

    for(Profile::vertex_descriptor_vector::const_iterator vit = mProfile.link().begin(); vit != mProfile.link().end(); ++ vit)
      this->add_vertex(B, *vit);

    for(Profile::Triangle_vector::const_iterator tit = mProfile.triangles().begin(); tit != mProfile.triangles().end(); ++ tit)
      this->add_triangle(B,*tit);

    B.end_surface();
  }
};

// Constructs a flat polyhedron containing just the link of a vertex
class Vertex_link_builder : public Link_builder
{
  const Surface& m_tm;
  Vertex_handle mV;

  const Surface& tm() { return m_tm; }

public:

  Vertex_link_builder(const Surface& aTM, Vertex_handle aV) : m_tm(aTM), mV(aV) {}

  void operator()(Surface::HalfedgeDS& hds)
  {
    Builder B(hds, true);

    B.begin_surface(1 + out_degree(mV,tm()), out_degree(mV,tm()));

    this->add_vertex(B,mV);
    Profile::Triangle_vector triangles;

    out_edge_iterator eb, ee;
    for(boost::tie(eb,ee) = halfedges_around_source(opposite(halfedge(mV,tm()),tm()),tm()); eb != ee; ++ eb)
    {
      halfedge_descriptor out_edge1 = *eb;
      halfedge_descriptor out_edge2 = out_edge1->opposite()->next();
      vertex_descriptor v1 = target(out_edge1,tm());
      vertex_descriptor v2 = target(out_edge2,tm());

      this->add_vertex(B,v1);

      triangles.push_back(Profile::Triangle(mV,v1,v2));
    }

    for(Profile::Triangle_vector::const_iterator tit = triangles.begin(); tit != triangles.end(); ++ tit)
      this->add_triangle(B,*tit);

    B.end_surface();
  }
};

SurfaceSP create_edge_link (const Profile& aProfile)
{
  SurfaceSP rLink(new Surface);

  try
  {
    Edge_link_builder lBuilder(aProfile);

    rLink->delegate(lBuilder);
  }
  catch(...)
  {
    cerr << "Unable to create local mesh" << endl;
  }

  return rLink;
}

SurfaceSP create_vertex_link (const Profile& aProfile, Vertex_handle aV)
{
  SurfaceSP rLink(new Surface);

  try
  {
    Vertex_link_builder lBuilder(aProfile.surface(),aV);

    rLink->delegate(lBuilder);
  }
  catch (...)
  {
    cerr << "Unable to create local mesh" << endl;
  }

  return rLink;
}

void write (SurfaceSP aSurface, string aName)
{
  ofstream out(aName.c_str());
  out << *aSurface;
}

template<class T>
string opt2str (const boost::optional<T>& o)
{
  ostringstream ss;
  if(o)
    ss << *o;
  else
    ss << "<none>";

  return ss.str();
}

template<class N>
string float2str (N n)
{
  return boost::str(boost::format("%+05.19g") % to_double(n));
}

template<class P>
string point2str (const P& p)
{
  return boost::str(boost::format("(%+05.19g %+05.19g %+05.19g)") % to_double(p.x()) % to_double(p.y()) % to_double(p.z()));
}

template<class P>
string optpoint2str (const boost::optional<P>& p)
{
  ostringstream ss;
  if(p)
       ss << point2str(*p);
  else ss << "<none>";

  return ss.str();
}
template<class N>
string optfloat2str (const boost::optional<N>& n)
{
  ostringstream ss;
  if(n)
       ss << float2str(*n);
  else ss << "<none>";

  return ss.str();
}

template<class V>
string vertex2str (const V& v)
{
  ostringstream ss;
  ss << "[V" << v->id() << point2str(v->point()) << "]";
  return ss.str();
}

template<class E>
string edge2str (const E& e)
{
  ostringstream ss;
  ss << "{E" << e->id() << vertex2str(e->opposite()->vertex()) << "->" << vertex2str(e->vertex()) << "}";
  return ss.str();
}

template<class T> ostream&  operator << (ostream& os, const boost::optional<T>& o) { return os << opt2str(o); }

string normalize_EOL (string line)
{
  string::size_type l = line.length();
  string::size_type d = (l > 0 && line[l-1] == '\r') ? 1 : 0;
  return line.substr(0, l-d);
}

#define REPORT_ERROR(msg) error(__FILE__,__LINE__,0, msg);

#define REPORT_ERROR2(pred,msg) REPORT_ERROR(msg)

#define CHECK_MSG(pred,msg) if(!(pred)) REPORT_ERROR2(#pred,msg)

#define CHECK(pred) CHECK_MSG(pred,string(""))

#define CHECK_EQUAL(x,y)       CHECK_MSG(((x)==(y)), str(format("Assertion failed: %1%(=%2%)==%3%(=%4%)") % (#x) % (x) % (#y) % (y)))
#define CHECK_NOT_EQUAL(x,y)   CHECK_MSG(((x)!=(y)), str(format("Assertion failed: %1%(=%2%)!=%3%(=%4%)") % (#x) % (x) % (#y) % (y)))

class Visitor : public SMS::Edge_collapse_visitor_base<Surface>
{

public :

  Visitor(string aTestCase) : mTestCase(aTestCase)
  {
#ifdef CGAL_SMS_TRACE_IMPL
    ::internal::cgal_enable_sms_trace = true;
#endif
    mStep = 0;
  }

  virtual ~Visitor()
  {
  }

  virtual void OnStarted(Surface&) const
  {
  }

  virtual void OnFinished (Surface& aSurface) const
  {
    CHECK(aSurface.is_valid());
  }

  virtual void OnCollected(const Profile& aProfile, const boost::optional<FT>& aCost) const
  {
    TEST_TRACE(str (format("Collecting %1% : cost=%2%") % edge2str(aProfile.v0_v1()) % optfloat2str(aCost)));
  }

  virtual void OnCollapsing(const Profile& aProfile, const boost::optional<Point>& aP) const
  {
    TEST_TRACE(str (format("S %1% - Collapsing %2% : placement=%3%") % mStep % edge2str(aProfile.v0_v1()) % optpoint2str(aP)));

    //mBefore = create_edge_link(aProfile);
  }

  virtual void OnCollapsed(const Profile& aProfile, const Vertex_handle& aV) const
  {
    // Some collapse can result in self-intersections in any case, so we can't mark it as a failure
    if(false && Is_self_intersecting(aProfile.surface()))
    {
      SurfaceSP lAfter = create_vertex_link(aProfile, aV);

      write(mBefore, str(format("%1%.step-%2%-before.off") % mTestCase % mStep));
      write(lAfter , str(format("%1%.step-%2%-after.off")  % mTestCase % mStep));

      REPORT_ERROR(str(format("Resulting surface self-intersects after step %1% (%2% edges left)") % mStep % (aProfile.surface().size_of_halfedges() / 2)));
    }

    ++mStep;
  }

private:

  void error (char const* file, int line, char const* pred, string msg) const
  {
    cerr << "ERROR in " << file << " at " << line << endl;
    if(pred)
      cerr << "  Assertion failed: " << pred << endl;
    cerr << "  " << msg << endl;

    throw runtime_error("");
  }

private:
  string            mTestCase;
  mutable int       mStep;
  mutable SurfaceSP mBefore;

};

bool sWriteResult = false;
bool sUseMP       = false;
int  sStop        = 1;

bool Test (string aName)
{
  bool rSucceeded = false;

  try
  {
    string off_name = aName;

    ifstream off_is(off_name.c_str());
    if(off_is)
    {
      Surface lSurface;
      off_is >> lSurface;
      if(lSurface.is_valid())
      {
        if(lSurface.is_pure_triangle())
        {
          cerr << "Processing " << aName << " (" << (lSurface.size_of_halfedges() / 2) << " edges)" << endl;

          Visitor vis(aName);

          set_halfedgeds_items_id(lSurface);

          SMS::Count_stop_predicate<Surface> stop(sStop);

          Real_timer t; t.start();

          if(sUseMP)
          {
            SMS::Edge_length_cost  <Surface> cost;
            SMS::Midpoint_placement<Surface> placement;

            edge_collapse(lSurface,stop,CGAL::parameters::get_cost(cost).get_placement(placement).visitor(vis));
          }
          else
          {
            SMS::LindstromTurk_cost     <Surface> cost;
            SMS::LindstromTurk_placement<Surface> placement;

            edge_collapse(lSurface,stop,CGAL::parameters::get_cost(cost).get_placement(placement).visitor(vis));
          }

          t.stop();

          rSucceeded = lSurface.is_valid() && lSurface.is_pure_triangle();

          cerr << "\r" << aName << (rSucceeded ? " succeeded" : "FAILED") << ". Elapsed time=" << t.time() << " seconds." << endl;

          if(sWriteResult)
          {
            string out_name = off_name;
            out_name.replace(off_name.find_last_of('.'), 15, ".simplified.off");
            ofstream out(out_name.c_str());
            out << lSurface;
            cerr << "Result written into: " << out_name << endl;
          }
        }
        else
        {
          cerr << "Surface is not triangulated (has faces with more than 3 sides): " << aName << endl;
        }
      }
      else
      {
        cerr << "Invalid surface: " << aName << endl;
      }
    }
    else
    {
      cerr << "Unable to open test file " << aName << endl;
    }
  }
  catch (const std::exception& x)
  {
    string what(x.what());
    if(what.length() > 0)
      cerr << "Exception caught: " << what << endl;
  }

  return rSucceeded;
}

int main(int argc, char** argv)
{
  set_error_handler  (error_handler);
  set_warning_handler(error_handler);

  cout << setprecision(4);
  cerr << setprecision(4);

  vector<string> lCases;

  for(int i = 1; i < argc; ++i)
  {

    if(argv[i][0] == '-')
    {
      switch(argv[i][1])
      {
        case 'w': sWriteResult = true; break;
        case 's': sStop = atoi(&argv[i][2]); break;
        case 'm': sUseMP = true; break;
      }
    }
    else
    {
      string c(normalize_EOL(string(argv[i])));
      string::size_type pos = c.find_last_of(".");
      if(pos != string::npos)
      {
        string ext = c.substr(pos);
        if(ext == ".off")
          lCases.push_back(c);
      }
    }
  }

  if(lCases.size() == 0)
  {
    cout << "collapse_edge_test file0 [-w] [-sNUM_EDGES] [-m] file1 ... fileN\n";
    return 1;
  }
  else
  {
    unsigned lOK = 0;
    for(vector<string>::const_iterator it = lCases.begin(); it != lCases.end(); ++ it)
    {
     if(Test(*it))
       ++ lOK;
    }

    cout << endl
         << lOK                    << " cases succedded." << endl
         << (lCases.size() - lOK) << " cases failed." << endl;

    return lOK == lCases.size() ? 0 : 1;
  }
}

#ifdef WHY_DO_THESE_DONT_EXIST_IN_CGAL_UNDER_LINUX
namespace CGAL
{

void assertion_fail(const char* expr, const char* file, int line)
{
  assertion_fail(expr,file,line,"");
}
void precondition_fail(const char* expr, const char* file, int line)
{
  precondition_fail(expr,file,line,"");
}
void postcondition_fail(const char* expr, const char* file, int line)
{
  postcondition_fail(expr,file,line,"");
}

}
#endif
