#ifndef LOCAL_TIME_TIME_ZONE_HPP__
#define LOCAL_TIME_TIME_ZONE_HPP__

/* Copyright (c) 2003-2004 CrystalClear Software, Inc.
 * Subject to the Boost Software License, Version 1.0. 
 * (See accompanying file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author: Jeff Garland, Bart Garst
 * $Date$
 */

#include "boost/date_time/time_zone_base.hpp"
#include "boost/date_time/time_zone_names.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/local_time/dst_transition_day_rules.hpp"
#include "boost/shared_ptr.hpp"

namespace boost {
namespace local_time {

  typedef boost::date_time::time_zone_names time_zone_names;
  typedef boost::date_time::dst_adjustment_offsets<boost::posix_time::time_duration> dst_adjustment_offsets;
  typedef boost::date_time::time_zone_base<boost::posix_time::ptime> time_zone_base;

  //! A real time zone
  class time_zone : public time_zone_base {
  public:
    typedef boost::posix_time::time_duration time_duration_type;
    typedef time_zone_base base_type;
    time_zone(const time_zone_names& zone_names,   
              const time_duration_type& base_utc_offset,
              const dst_adjustment_offsets& dst_offset,
              boost::shared_ptr<dst_calc_rule> calc_rule) :
      zone_names_(zone_names),
      has_dst_(calc_rule != 0),
      base_utc_offset_(base_utc_offset),
      dst_offsets_(dst_offset),
      dst_calc_rules_(calc_rule)
    {};
    virtual ~time_zone() {};
    virtual std::string dst_zone_abbrev() const
    {
      return zone_names_.dst_zone_abbrev();
    }
    virtual std::string std_zone_abbrev() const
    {
      return zone_names_.std_zone_abbrev();
    }
    virtual std::string dst_zone_name() const
    {
      return zone_names_.dst_zone_name();
    }
    virtual std::string std_zone_name() const
    {
      return zone_names_.std_zone_name();
    }
    //! True if zone uses daylight savings adjustments
    virtual bool has_dst() const
    {
      return has_dst_;
    }
    //! Local time that DST starts -- undefined if has_dst is false
    virtual posix_time::ptime dst_local_start_time(gregorian::greg_year y) const
    {
      gregorian::date d(1900,1,1);
      if (dst_calc_rules_) {
        d = dst_calc_rules_->start_day(y);
      }
      return posix_time::ptime(d, dst_offsets_.dst_start_offset_);
    }
    //! Local time that DST ends -- undefined if has_dst is false
    virtual posix_time::ptime dst_local_end_time(gregorian::greg_year y) const
    {
      gregorian::date d(1900,1,1);
      if (dst_calc_rules_) {
        d = dst_calc_rules_->end_day(y);
      }
      return posix_time::ptime(d, dst_offsets_.dst_end_offset_);
    }
    //! Base offset from UTC for zone (eg: -07:30:00)
    virtual time_duration_type base_utc_offset() const
    {
      return base_utc_offset_;
    }
    //! Adjustment forward or back made while DST is in effect
    virtual time_duration_type dst_offset() const
    {
      return dst_offsets_.dst_adjust_;
    }
  private:
    time_zone_names zone_names_;
    bool has_dst_;
    time_duration_type base_utc_offset_;
    dst_adjustment_offsets dst_offsets_;
    boost::shared_ptr<dst_calc_rule> dst_calc_rules_;
  };

} }//namespace



#endif
