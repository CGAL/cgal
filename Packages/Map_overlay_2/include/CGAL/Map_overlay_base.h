#ifndef CGAL_MAP_OVERLAY_BASE_H
#define CGAL_MAP_OVERLAY_BASE_H

CGAL_BEGIN_NAMESPACE

template <class Arrangement_, class Map_overlay_ChangeNotification_>
class Map_overlay_base
{
public:
  typedef Arrangement_                     Arrangement;
  typedef Map_overlay_ChangeNotification_  Map_overlay_ChangeNotification;

  Map_overlay_base() {}
  
  virtual void map_overlay(const Arrangement &a1, 
                           const Arrangement &a2, 
                           Map_overlay_ChangeNotification *pmwx_change_notf, 
                           Arrangement &result) = 0;

  virtual ~Map_overlay_base() {}
};

CGAL_END_NAMESPACE

#endif
