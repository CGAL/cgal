#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_STORAGE_TRAITS_WITH_INFO_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_STORAGE_TRAITS_WITH_INFO_2_H 1

#include <CGAL/basic.h>
#include <CGAL/Storage_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_site_with_info_2.h>

CGAL_BEGIN_NAMESPACE

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<class STraits>
class Construct_storage_site_with_info_2
  : public Construct_storage_site_2<STraits>
{
public:
  typedef STraits                                    Storage_traits;
  typedef typename Storage_traits::Storage_site_2    Storage_site_2;
  typedef typename Storage_traits::Point_handle      Point_handle;

  typedef Storage_site_2                             result_type;
  //  struct Arity {};

protected:
  typedef Construct_storage_site_2<Storage_traits>   Base;
  typedef typename Storage_traits::Info              Info;
  typedef typename Storage_traits::Convert_info      Convert_info;
  typedef typename Storage_traits::Merge_info        Merge_info;

public:
  // constructs the point of intersection
  inline
  result_type operator()(const Storage_site_2& ss0,
			 const Storage_site_2& ss1) const {
    Storage_site_2 ssx = Base::operator()(ss0, ss1);
    Info infox = Merge_info()(ss0.info(), ss1.info());
    ssx.set_info(infox);
    return ssx;
  }

  // constructs the subsegment with supporting segment ss0 and
  // endpoints the point of intersection of ss1 and ss0; the boolean
  // determines if the first or segment subsegment is constructed
  inline
  result_type operator()(const Storage_site_2& ss0,
			 const Storage_site_2& ss1,
			 bool first) const {
    Storage_site_2 s = Base::operator()(ss0, ss1, first);
    Info is = Convert_info()(ss0.info(), ss1.info(), first);
    s.set_info(is);
    return s;
  }

  using Base::operator();
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<typename Info_>
struct Dummy_convert_info
{
  typedef Info_  Info;
  typedef Info   result_type;

  inline
  Info operator()(const Info&, const Info&, bool) const {
    return Info();
  }
};

template<class Info_>
struct Dummy_merge_info
{
  typedef Info_  Info;
  typedef Info   result_type;

  inline
  Info operator()(const Info&, const Info&) const {
    return Info();
  }
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<class STraits_base,
	 typename Info_ = char,
	 class Converter = Dummy_convert_info<Info_>,
	 class Merger = Dummy_convert_info<Info_> >
class Storage_traits_with_info_2
  : public STraits_base
{
public:
  typedef Info_                                    Info;
  typedef Converter                                Convert_info;
  typedef Merger                                   Merge_info;

private:
  typedef STraits_base                             Base;
  typedef typename Base::Storage_site_2            Base_storage_site_2;

  typedef Storage_traits_with_info_2<Base,Info,Convert_info,Merge_info> Self;

public:
  typedef typename Base::Geom_traits               Geom_traits;
#if 0
  typedef typename Base::Point_2                   Point_2;
  typedef typename Base::Site_2                    Site_2;
  typedef typename Base::Point_container           Point_container;
  typedef typename Base::Point_handle              Point_handle;
#endif

  typedef
  Segment_Delaunay_graph_storage_site_with_info_2<Self,Info,Base_storage_site_2>
  Storage_site_2;

  typedef Construct_storage_site_with_info_2<Self>
  Construct_storage_site_2;

  // MK::FIGURE OUT HOW TO PASS A REFERENCE TO GEOM_TRAITS AND HAVE
  // DEFAULT CONSTRUCTOR AS WELL IF POSSIBLE
  Storage_traits_with_info_2(const Geom_traits& gt = Geom_traits())
    : Base(gt) {}

  inline Construct_storage_site_2
  construct_storage_site_2_object() const {
    return Construct_storage_site_2();
  }

  inline Convert_info
  convert_info_object() const {
    return Convert_info();
  }

  inline Merge_info
  merge_info_object() const {
    return Merge_info();
  }
};



CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_STORAGE_TRAITS_WITH_INFO_2_H
