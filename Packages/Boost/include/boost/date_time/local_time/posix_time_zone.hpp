#ifndef _DATE_TIME_POSIX_TIME_ZONE__
#define _DATE_TIME_POSIX_TIME_ZONE__ 

/* Copyright (c) 2003-2004 CrystalClear Software, Inc.
 * Subject to the Boost Software License, Version 1.0. (See accompanying
 * file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author: Jeff Garland, Bart Garst
 * $Date$
 */

#include <string>
#include <sstream>
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/time_zone_names.hpp"
#include "boost/date_time/time_zone_base.hpp"
#include "boost/date_time/local_time/dst_transition_day_rules.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/tokenizer.hpp"
#include <stdexcept>

namespace boost{
namespace local_time{

  //! simple exception for UTC and Daylight savings start/end offsets
  struct bad_offset : public std::out_of_range
  {
    bad_offset(std::string _msg="") : std::out_of_range(std::string("Offset out of range: " + _msg)) {}
  };
  //! simple exception for UTC daylight savings adjustment
  struct bad_adjustment : public std::out_of_range
  {
    bad_adjustment(std::string _msg="") : std::out_of_range(std::string("Adjustment out of range: " + _msg)) {}
  };
  
  typedef boost::date_time::time_zone_names time_zone_names;
  typedef boost::date_time::dst_adjustment_offsets<boost::posix_time::time_duration> dst_adjustment_offsets;
  typedef boost::date_time::time_zone_base<boost::posix_time::ptime> time_zone_base;

  //! A time zone class constructed from a POSIX time zone string
  /*! A POSIX time zone string takes to form of:<br>
   * "std offset dst [offset],start[/time],end[/time]" (w/no spaces)
   * 'std' specifies the abbrev of the time zone.<br> 
   * 'offset' is the offset from UTC.<br>
   * 'dst' specifies the abbrev of the time zone during daylight savings time.<br>
   * The second offset is how many hours changed during DST. Default=1<br>
   * 'start' & 'end' are the dates when DST goes into (and out of) effect.<br>
   * 'offset' takes the form of: [+|-]hh[:mm[:ss]] {h=0-23, m/s=0-59}<br>
   * 'time' and 'offset' take the same form. Time defaults=02:00:00<br>
   * 'start' & 'end' can be one of three forms:<br>
   * Mm.w.d {month=1-12, week=1-5 (5 is always last), day=0-6}<br>
   * Jn {n=1-365 Feb29 is never counted}<br>
   * n  {n=0-365 Feb29 is counted in leap years}<br>
   * Example "PST-5PDT01:00:00,M4.1.0/02:00:00,M10.1.0/02:00:00"
   * <br>
   * Exceptions will be thrown under these conditions:<br>
   * An invalid date spec (see date class)<br>
   * A boost::local_time::bad_offset exception will be thrown for:<br>
   * A DST start or end offset that is negative or more than 24 hours<br>
   * A UTC zone that is greater than +12 or less than -12 hours<br>
   * A boost::local_time::bad_adjustment exception will be thrown for:<br>
   * A DST adjustment that is 24 hours or more (positive or negative)<br>
   */
  class posix_time_zone : public time_zone_base  {
  public:
    typedef boost::posix_time::time_duration time_duration_type;
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

    //! Construct from a POSIX time zone string
    posix_time_zone(const std::string& s) : 
      zone_names_("std_name","std_abbrev","no-dst","no-dst"),
      has_dst_(false), 
      base_utc_offset_(posix_time::hours(0)),
      dst_offsets_(posix_time::hours(0),posix_time::hours(0),posix_time::hours(0)),
      dst_calc_rules_()
    {
      boost::char_separator<char> sep(",");
      tokenizer tokens(s, sep);
      tokenizer::iterator it = tokens.begin();
      calc_zone(*it++);
      if(has_dst_){
        std::string tmp_str = *it++;
        calc_rules(tmp_str, *it);
      }
    } 
    virtual ~posix_time_zone() {};
    //!String for the zone when not in daylight savings (eg: EST)
    virtual std::string std_zone_abbrev()const
    {
      return zone_names_.std_zone_abbrev();
    }
    //!String for the timezone when in daylight savings (eg: EDT)
    /*! For those time zones that have no DST, an empty string is used */
    virtual std::string dst_zone_abbrev() const
    {
      return zone_names_.dst_zone_abbrev();
    }
    //!String for the zone when not in daylight savings (eg: Eastern Standard Time)
    /*! The full STD name is not extracted from the posix time zone string. 
     * Therefore, the STD abbreviation is used in it's place */
    virtual std::string std_zone_name()const
    {
      return zone_names_.std_zone_name();
    }
    //!String for the timezone when in daylight savings (eg: Eastern Daylight Time)
    /*! The full DST name is not extracted from the posix time zone string. 
     * Therefore, the STD abbreviation is used in it's place. For time zones 
     * that have no DST, an empty string is used */
    virtual std::string dst_zone_name()const
    {
      return zone_names_.dst_zone_name();
    }
    //! True if zone uses daylight savings adjustments otherwise false
    virtual bool has_dst()const
    {
      return has_dst_;
    }
    //! Local time that DST starts -- undefined if has_dst is false
    virtual posix_time::ptime dst_local_start_time(gregorian::greg_year y)const
    {
      gregorian::date d(1900,1,1);
      if(has_dst_)
      {
        d = dst_calc_rules_->start_day(y);
      }
      return posix_time::ptime(d, dst_offsets_.dst_start_offset_);
    }
    //! Local time that DST ends -- undefined if has_dst is false
    virtual posix_time::ptime dst_local_end_time(gregorian::greg_year y)const
    {
      gregorian::date d(1900,1,1);
      if(has_dst_)
      {
        d = dst_calc_rules_->end_day(y);
      }
      return posix_time::ptime(d, dst_offsets_.dst_end_offset_);
    }
    //! Base offset from UTC for zone (eg: -07:30:00)
    virtual time_duration_type base_utc_offset()const
    {
      return base_utc_offset_;
    }
    //! Adjustment forward or back made while DST is in effect
    virtual time_duration_type dst_offset()const
    {
      return dst_offsets_.dst_adjust_;
    }

  private:
    time_zone_names zone_names_;
    bool has_dst_;
    time_duration_type base_utc_offset_;
    dst_adjustment_offsets dst_offsets_;
    boost::shared_ptr<dst_calc_rule> dst_calc_rules_;

    /*! Extract time zone abbreviations for STD & DST as well
     * as the offsets for the time the shift occurs and how
     * much of a shift. At this time full time zone names are
     * NOT extracted so the abbreviations are used in their place */
    void calc_zone(const std::string& obj){
      std::stringstream ss("");
      std::string::const_iterator sit = obj.begin();
      std::string std_zone_abbrev("std_abbrev"), dst_zone_abbrev("");

      // get 'std' name/abbrev
      while(isalpha(*sit)){
        ss << *sit++;
      }
      std_zone_abbrev = ss.str(); 
      ss.str("");

      // get UTC offset
      if(sit != obj.end()){
        // get duration
        while(!isalpha(*sit) && sit != obj.end()){
        ss << *sit++;
    }
    base_utc_offset_ = posix_time::duration_from_string(ss.str()); 
    ss.str("");

    // base offset must be within range of -12 hours to +12 hours
    if(base_utc_offset_ < time_duration_type(-12,0,0) ||
        base_utc_offset_ > time_duration_type(12,0,0))
    {
        throw bad_offset(posix_time::to_simple_string(base_utc_offset_));
    }
      }

      // get DST data if given
      if(sit != obj.end()){
    has_dst_ = true;
    
        // get 'dst' name/abbrev
        while(isalpha(*sit)){
          ss << *sit++;
        }
        dst_zone_abbrev = ss.str(); 
        ss.str("");

        // get DST offset if given
        if(sit != obj.end()){
          // get duration
          while(!isalpha(*sit) && sit != obj.end()){
            ss << *sit++;
          }
          dst_offsets_.dst_adjust_ = 
                posix_time::duration_from_string(ss.str());
        ss.str("");
        }
        else{ // default DST offset
          dst_offsets_.dst_adjust_ = posix_time::hours(1);
        }

    // adjustment must be within +|- 1 day
    if(dst_offsets_.dst_adjust_ <= time_duration_type(-24,0,0) ||
        dst_offsets_.dst_adjust_ >= time_duration_type(24,0,0))
    {
      throw bad_adjustment(posix_time::to_simple_string(dst_offsets_.dst_adjust_));
    }
      }
      // full names not extracted so abbrevs used in their place
      zone_names_ = time_zone_names(std_zone_abbrev, std_zone_abbrev, dst_zone_abbrev, dst_zone_abbrev);
    }

    void calc_rules(const std::string& start, const std::string& end){
      boost::char_separator<char> sep("/");
      tokenizer st_tok(start, sep);
      tokenizer et_tok(end, sep);
      tokenizer::iterator sit = st_tok.begin();
      tokenizer::iterator eit = et_tok.begin();

      // generate date spec
      char x = std::string(*sit).at(0);
      if(x == 'M'){
        M_func(*sit, *eit);
      }
      else if(x == 'J'){
        julian_no_leap(*sit, *eit);
      }
      else{
        julian_day(*sit, *eit);
      }

      ++sit;
      ++eit;
      // generate durations
      // starting offset
      if(sit != st_tok.end()){
        dst_offsets_.dst_start_offset_ = posix_time::duration_from_string(*sit);
      }
      else{
        // default
        dst_offsets_.dst_start_offset_ = posix_time::hours(2);
      }
      // start/end offsets must fall on given date
      if(dst_offsets_.dst_start_offset_ < time_duration_type(0,0,0) ||
          dst_offsets_.dst_start_offset_ >= time_duration_type(24,0,0))
      {
        throw bad_offset(posix_time::to_simple_string(dst_offsets_.dst_start_offset_));
      }

      // ending offset
      if(eit != et_tok.end()){
        dst_offsets_.dst_end_offset_ = posix_time::duration_from_string(*eit);
      }
      else{
        // default
        dst_offsets_.dst_end_offset_ = posix_time::hours(2);
      }
      // start/end offsets must fall on given date
      if(dst_offsets_.dst_end_offset_ < time_duration_type(0,0,0) ||
        dst_offsets_.dst_end_offset_ >= time_duration_type(24,0,0))
      {
        throw bad_offset(posix_time::to_simple_string(dst_offsets_.dst_end_offset_));
      }
    }

    /* Parses out a start/end date spec from a posix time zone string.
     * Date specs come in three possible formats, this function handles
     * the 'M' spec. Ex "M2.2.4" => 2nd month, 2nd week, 4th day .
     */
    void M_func(const std::string& s, const std::string& e){
      typedef gregorian::nth_kday_of_month nkday;
      int sm=0,sw=0,sd=0,em=0,ew=0,ed=0; // start/end month,week,day
      char_separator<char> sep("M.");
      tokenizer stok(s, sep), etok(e, sep);
      
      tokenizer::iterator it = stok.begin();
      sm = lexical_cast<int>(*it++);
      sw = lexical_cast<int>(*it++);
      sd = lexical_cast<int>(*it);
     
      it = etok.begin();
      em = lexical_cast<int>(*it++);
      ew = lexical_cast<int>(*it++);
      ed = lexical_cast<int>(*it);

      dst_calc_rules_ = shared_ptr<dst_calc_rule>(
        new nth_kday_dst_rule(
          nth_last_dst_rule::start_rule(
            static_cast<nkday::week_num>(sw),sd,sm), 
          nth_last_dst_rule::start_rule(
            static_cast<nkday::week_num>(ew),ed,em) 
          )
      );
    }
    
    //! Julian day. Feb29 is never counted, even in leap years
    // expects range of 1-365
    void julian_no_leap(const std::string& s, const std::string& e){
      typedef gregorian::gregorian_calendar calendar;
      const unsigned short year = 2001; // Non-leap year
      int sm=1, sd=0;
      sd = lexical_cast<int>(s.substr(1)); // skip 'J'
      while(sd >= calendar::end_of_month_day(year,sm)){
        sd -= calendar::end_of_month_day(year,sm++);
      }
      int em=1, ed=0;
      ed = lexical_cast<int>(e.substr(1)); // skip 'J'
      while(ed >= calendar::end_of_month_day(year,em)){
        ed -= calendar::end_of_month_day(year,em++);
      }

      dst_calc_rules_ = shared_ptr<dst_calc_rule>(
        new partial_date_dst_rule(
          partial_date_dst_rule::start_rule(
            sd, static_cast<date_time::months_of_year>(sm)), 
          partial_date_dst_rule::end_rule(
            ed, static_cast<date_time::months_of_year>(em)) 
          )
      );
    }

    //! Julian day. Feb29 is always counted, but exception thrown in non-leap years
    // expects range of 0-365
    void julian_day(const std::string& s, const std::string& e){
      int sd=0, ed=0;
      sd = lexical_cast<int>(s);
      ed = lexical_cast<int>(e);
      dst_calc_rules_ = shared_ptr<dst_calc_rule>(
        new partial_date_dst_rule(
          partial_date_dst_rule::start_rule(++sd),// args are 0-365
          partial_date_dst_rule::end_rule(++ed) // pd expects 1-366
          )
      );
    }

  };

} } // namespace boost::local_time


#endif // _DATE_TIME_POSIX_TIME_ZONE__ 
