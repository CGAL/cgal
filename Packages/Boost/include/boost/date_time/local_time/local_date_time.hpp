#ifndef LOCAL_TIME_LOCAL_DATE_TIME_HPP__
#define LOCAL_TIME_LOCAL_DATE_TIME_HPP__

/* Copyright (c) 2003-2004 CrystalClear Software, Inc.
 * Subject to the Boost Software License, Version 1.0. 
 * (See accompanying file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author: Jeff Garland, Bart Garst
 * $Date$
 */


#include "boost/date_time/time.hpp"
#include "boost/date_time/posix_time/posix_time.hpp" //todo remove?
#include "boost/shared_ptr.hpp"
#include "boost/date_time/dst_rules.hpp"
#include "boost/date_time/time_zone_base.hpp"
#include "boost/date_time/special_defs.hpp"
#include <sstream>

namespace boost {
namespace local_time {

  //! simple exception for reporting when STD or DST cannot be determined
  struct ambiguous_result : public std::logic_error
  {
    ambiguous_result (std::string _msg="") :
      std::logic_error(std::string("Daylight Savings Results are ambiguous: " + _msg)) {}
  };
  // 
  struct InvalidTimeLabel : public std::logic_error
  {
    InvalidTimeLabel(std::string _msg="") :
      std::logic_error(std::string("Time label given is invalid: " + _msg)) {}
  };
  struct DSTNotValid: public std::logic_error
  {
    DSTNotValid(std::string _msg="") :
      std::logic_error(std::string("is_dst flag does not match resulting dst for time label given: " + _msg)) {}
  };

  //TODO: I think these should be in local_date_time_base and not 
  // necessarily brought into the namespace
  using date_time::time_is_dst_result;
  using date_time::is_in_dst;
  using date_time::is_not_in_dst;
  using date_time::ambiguous;
  using date_time::invalid_time_label;
 
  template<class utc_time_=posix_time::ptime, 
    class tz_type=date_time::time_zone_base<typename posix_time::ptime> >
  class local_date_time_base :  public date_time::base_time<boost::posix_time::ptime, 
                                                            boost::posix_time::posix_time_system> { 
  public:
    typedef utc_time_ utc_time_type;
    typedef typename utc_time_type::time_duration_type time_duration_type;
    typedef typename utc_time_type::date_type date_type;
    typedef typename date_type::duration_type date_duration_type;
    typedef typename utc_time_type::time_system_type time_system_type;
    /** This constructor interprets the passed time as a UTC time.
     *  So, for example, if the passed timezone is UTC-5 then the
     *  time will be adjusted back 5 hours.  The time zone allows for 
     *  automatic calculation of whether the particular time is adjusted for
     *  daylight savings, etc.
     *  If the time zone shared pointer is null then time stays unadjusted.
     *@param t A UTC time
     *@param tz Timezone for to adjust the UTC time to.
     */
    local_date_time_base(utc_time_type t, 
                         boost::shared_ptr<tz_type> tz) : 
      date_time::base_time<utc_time_type, time_system_type>(t),
      zone_(tz)
    {
      // param was already utc so nothing more to do
    } 

    /** This constructs a local time -- the passed time information 
     * understood to be in the passed tz. The DST flag must be passed 
     * to indicate whether the time is in daylight savings or not.  
     *  @throws -- InvalidTimeLabel if the time passed does not exist in 
     *             the given locale. The non-existant case occurs typically 
     *             during the shift-back from daylight savings time.  When 
     *             the clock is shifted forward a range of times 
     *             (2 am to 3 am in the US) is skipped and hence is invalid.
     *  @throws -- DSTNotValid if the DST flag is passed for a period 
     *             where DST is not active.
     */
    local_date_time_base(date_type d, 
                         time_duration_type td,
                         boost::shared_ptr<tz_type> tz, 
                         bool dst_flag) : //necessary for constr_adj()
      date_time::base_time<utc_time_type,time_system_type>(construction_adjustment(utc_time_type(d, td), tz, dst_flag)),
      zone_(tz)
    {
      if(tz != NULL && tz->has_dst()){
        
        // d & td are already local so we use them
        time_is_dst_result result = check_dst(d, td, tz);
        bool in_dst = (result == is_in_dst); // less processing than is_dst()

        // ambig occurs at end, invalid at start
        if(result == invalid_time_label){
          // Ex: 2:15am local on trans-in day in nyc, dst_flag irrelevant
          std:: stringstream ss;
          ss << "time given: " << d << ' ' << td;
          throw InvalidTimeLabel(ss.str());
        }
        else if(result != ambiguous && in_dst != dst_flag){
          // is dst_flag accurate?
          // Ex: false flag in NYC in June
          std:: stringstream ss;
          ss << "flag given: " << (dst_flag ? "dst=true" : "dst=false")
            << ", dst calculated: " << (in_dst ? "dst=true" : "dst=false");
          throw DSTNotValid(ss.str());
        }
        
        // everything checks out and conversion to utc already done
      }
    } 
    
    //TODO maybe not the right set...Ignore the last 2 for now...
    enum DST_CALC_OPTIONS { EXCEPTION_ON_ERROR, NOT_DATE_TIME_ON_ERROR }; 
                            //ASSUME_DST_ON_ERROR, ASSUME_NOT_DST_ON_ERROR };

    /** This constructs a local time -- the passed time information 
     * understood to be in the passed tz.  The DST flag is calculated 
     * according to the specified rule. 
     */
    local_date_time_base(date_type d, 
                         time_duration_type td,
                         boost::shared_ptr<tz_type> tz, 
                         DST_CALC_OPTIONS calc_option) : 
      // dummy value - time_ is set in constructor code
      date_time::base_time<utc_time_type,time_system_type>(utc_time_type(d,td)),
      zone_(tz)
    {
      time_is_dst_result result = check_dst(d, td, tz);
      if(result == ambiguous) {
        if(calc_option == EXCEPTION_ON_ERROR){
          std:: stringstream ss;
          ss << "time given: " << d << ' ' << td;
          throw ambiguous_result(ss.str());
        }
        else{ // NADT on error
          time_ = posix_time::posix_time_system::get_time_rep(date_type(date_time::not_a_date_time), time_duration_type(date_time::not_a_date_time));
        }
      }
      else if(result == invalid_time_label){
        if(calc_option == EXCEPTION_ON_ERROR){
          std:: stringstream ss;
          ss << "time given: " << d << ' ' << td;
          throw InvalidTimeLabel(ss.str());
        }
        else{ // NADT on error
          time_ = posix_time::posix_time_system::get_time_rep(date_type(date_time::not_a_date_time), time_duration_type(date_time::not_a_date_time));
        }
      }
      else if(result == is_in_dst){
        utc_time_type t = 
          construction_adjustment(utc_time_type(d, td), tz, true);
        time_ = posix_time::posix_time_system::get_time_rep(t.date(), 
                                                            t.time_of_day());
      }
      else{
        utc_time_type t = 
          construction_adjustment(utc_time_type(d, td), tz, false);
        time_ = posix_time::posix_time_system::get_time_rep(t.date(), 
                                                            t.time_of_day());
      }
    }


    //! takes a date and time_duration representing a local time
    /*! 
     * returns TODO: ... 
     */
    static time_is_dst_result check_dst(date_type d, 
                                        time_duration_type td,
                                        boost::shared_ptr<tz_type> tz) 
    {
      if(tz != NULL && tz->has_dst()) {
        typedef typename date_time::dst_calculator<date_type, time_duration_type> dst_calculator;
        return dst_calculator::local_is_dst(
            d, td, 
            tz->dst_local_start_time(d.year()).date(), 
            tz->dst_local_start_time(d.year()).time_of_day(), 
            tz->dst_local_end_time(d.year()).date(), 
            tz->dst_local_end_time(d.year()).time_of_day(),
            tz->dst_offset()
        );
      }
      else{
        return is_not_in_dst;
      }
    }

    //! Simple destructor, releases time zone if last referer
    ~local_date_time_base() {};

    //! Copy constructor
    local_date_time_base(const local_date_time_base& rhs) : 
      date_time::base_time<utc_time_type, time_system_type>(rhs),
      zone_(rhs.zone_)
    {}

    //! Special values constructor
    explicit local_date_time_base(const boost::date_time::special_values sv,
                                  boost::shared_ptr<tz_type> tz = boost::shared_ptr<tz_type>()) :
      date_time::base_time<utc_time_type, time_system_type>(utc_time_type(sv)),
      zone_(tz)
    {}

    boost::shared_ptr<tz_type> zone() const 
    {
      return zone_;
    }
    //! returns false is time_zone is NULL and if time value is a special_value
    bool is_dst() const
    {
      if(zone_ != NULL && zone_->has_dst() && !is_special()) {
        // check_dst takes a local time, *this is utc
        utc_time_type lt(time_);
        lt += zone_->base_utc_offset();
        // dst_offset only needs to be considered with ambiguous time labels
        // make that adjustment there

        switch(check_dst(lt.date(), lt.time_of_day(), zone_)){
          case is_not_in_dst:
            return false;
            break;
          case is_in_dst:
            return true;
            break;
          case ambiguous: 
            if(lt + zone_->dst_offset() < zone_->dst_local_end_time(lt.date().year())) {
              return true;
            }
            break;
          case invalid_time_label: 
            if(lt >= zone_->dst_local_start_time(lt.date().year())) {
              return true;
            }
            break;
        }
      }
      return false;
    }
    //! Returns object's time value as a utc representation
    utc_time_type utc_time() const 
    {
      return utc_time_type(time_);
    }
    //! Returns object's time value as a local representation
    utc_time_type local_time() const 
    {
      if(zone_ != NULL){
        utc_time_type lt = this->utc_time() + zone_->base_utc_offset();
        if (is_dst()) {
          lt += zone_->dst_offset();
        }
        return lt;
      }
      return utc_time_type(time_);
    }
    //! Returns string in the form "2003-Aug-20 05:00:00 EDT"
    /*! Returns string in the form "2003-Aug-20 05:00:00 EDT". If
     * time_zone is NULL the time zone abbreviation will be "UTC". The time 
     * zone abbrev will not be included if calling object is a special_value*/
    std::string to_string() const
    {
      std::stringstream ss;
      if(is_special()){
        ss << utc_time();
        return ss.str();
      }
      if(zone_ == NULL) {
        ss << utc_time() << " UTC";
        return ss.str();
      }
      bool is_dst_ = is_dst();
      utc_time_type lt = this->utc_time() + zone_->base_utc_offset();
      if (is_dst_) {
        lt += zone_->dst_offset();
      }
      if(zone_ == NULL) {
        ss << utc_time() << " UTC";
        return ss.str();
      }
      ss << local_time() << " ";
      if (is_dst()) {
        ss << zone_->dst_zone_abbrev();
      }
      else {
        ss << zone_->std_zone_abbrev();
      }
      return ss.str();
    }
    /*! returns a local_date_time_base in the given time zone with the 
     * optional time_duration added. */
    local_date_time_base local_time_in(boost::shared_ptr<tz_type> new_tz, 
                                       time_duration_type td=time_duration_type(0,0,0)) const 
    {
      return local_date_time_base(utc_time_type(time_) + td, new_tz);
    }
 
      
    //! Equality comparison operator
    /*bool operator==(const date_time::base_time<boost::posix_time::ptime,boost::posix_time::posix_time_system>& rhs) const
    { // fails due to rhs.time_ being protected
      return date_time::base_time<boost::posix_time::ptime,boost::posix_time::posix_time_system>::operator==(rhs);
      //return this->time_ == rhs.time_;
    }*/
    //! Equality comparison operator
    bool operator==(const local_date_time_base& rhs) const
    {
      return time_system_type::is_equal(this->time_, rhs.time_);
    }
    //! Non-Equality comparison operator
    bool operator!=(const local_date_time_base& rhs) const
    {
      return !(*this == rhs);
    }
    //! Less than comparison operator
    bool operator<(const local_date_time_base& rhs) const
    {
      return time_system_type::is_less(this->time_, rhs.time_);
    }
    //! Less than or equal to comparison operator
    bool operator<=(const local_date_time_base& rhs) const
    {
      return (*this < rhs || *this == rhs);
    }
    //! Greater than comparison operator
    bool operator>(const local_date_time_base& rhs) const
    {
      return !(*this <= rhs);
    }
    //! Greater than or equal to comparison operator
    bool operator>=(const local_date_time_base& rhs) const
    {
      return (*this > rhs || *this == rhs);
    }

    /* Math operators. local_date_time_base inherits it's math operators 
     * from base_time. The shortcut operators (-= & +=) work but the 
     * others don't (they return a time_type which in local_date_time 
     * is a ptime). */
    //! Local_date_time + date_duration 
    local_date_time_base operator+(const date_duration_type& dd) const
    {
      return local_date_time_base(time_system_type::add_days(time_,dd), zone_);
    }
    //! Local_date_time - date_duration 
    local_date_time_base operator-(const date_duration_type& dd) const
    {
      return local_date_time_base(time_system_type::subtract_days(time_,dd), zone_);
    }
    //! Local_date_time + time_duration 
    local_date_time_base operator+(const time_duration_type& td) const
    {
      return local_date_time_base(time_system_type::add_time_duration(time_,td), zone_);
    }
    //! Local_date_time - time_duration 
    local_date_time_base operator-(const time_duration_type& td) const
    {
      return local_date_time_base(time_system_type::subtract_time_duration(time_,td), zone_);
    }
  private:
    boost::shared_ptr<tz_type> zone_;
    //bool is_dst_;

    /*! Adjust the passed in time to UTC?
     */
    utc_time_type construction_adjustment(utc_time_type t, 
                                          boost::shared_ptr<tz_type> zone,
                                          bool is_dst)
    {
      if(zone != NULL) {
        if(is_dst && zone->has_dst()) {
          t -= zone->dst_offset();
        } // else no adjust
        t -= zone->base_utc_offset();
      }
      return t;
    }
  };

  //!Use the default parameters to define local_date_time
  typedef local_date_time_base<> local_date_time;

} }


#endif
