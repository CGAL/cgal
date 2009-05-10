#include <vector>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

class Visitor
{
private :
  
  struct Data
  {
    Data ( size_t i, optional<NT> const& c, optional<Point> const& p )
      : 
       id(i)
      ,cost(c) 
      ,placement(p)
    {}
    
    size_t           id ;
    optional<NT>     cost ;
    optional<Point>  placement ;
  } ;
  typedef boost::shared_ptr<Data> Data_ptr ;
  
  typedef map<int,Data_ptr> Table ;
 
  
private :
  
  typedef char_separator<char> Separator ;
  typedef tokenizer<Separator> Tokenizer ;
  

  void ReadAudit( istream& in )
  {
    size_t section = 0 ;
    string line ;
    while ( getline(in,line) )
    {
      line = normalize_EOL(line);
      if ( line.length() > 0 )
      {
        TRACE( str ( format("AUDIT: %1%") % line ) ) ;
       
        switch ( section )
        {
          case 0 : epsilon_cost   = compute_epsilon_from_smallest_sample(toNT(line)); break ;
          case 1 : epsilon_sqdist = compute_epsilon_from_smallest_sample(toNT(line)); break ;
          
          default : 
     
            Tokenizer tk(line,Separator(" "));
            vector<string> tokens (tk.begin(),tk.end());

            CHECK_MSG( tokens.size() >= 1, str(format("Invalid audit line of type I, id field missing: %1%") % line) ) ;
            CHECK_MSG( tokens.size() >= 2, str(format("Invalid audit line of type I, cost field missing: %1%") % line) ) ;
            CHECK_MSG( tokens.size() >= 5, str(format("Invalid audit line of type I, placement field missing: %1%") % line) ) ;
            
            size_t id = lexical_cast<size_t>(tokens[0]);

            CHECK_MSG(audit_table.find(id)==audit_table.end()
                     ,str(format("Invalid audit line of type I, id field with duplicate value: %1%") % line)
                     );

                           
            optional<NT> cost ;
                  
            if ( tokens[1][0] != '?' )
              cost = toNT(tokens[1]);
                  
            optional<Point> placement ;
            if ( tokens.size() > 3 )
            {
              if ( tokens[2][0] != '?' )
              {
                NT x = toNT(tokens[2]);
                NT y = toNT(tokens[3]);
                NT z = toNT(tokens[4]);
                      
                placement = Point(x,y,z);
              }              
            }
            
            audit_table.insert(make_pair(id,Data_ptr( new Data(id, cost, placement) ) ) ) ;
            
            break ;
        }
        
        ++ section ;
      }  
    }
  }
  
public :
  
  Visitor ( string audit_name ) : infinite_cost(1e+8)
  {
    TRACE( str ( format("AUDIT FILE: %1%") % audit_name ) ) ;
    ifstream in(audit_name.c_str());  
    if ( in )
         ReadAudit(in);
    else REPORT_ERROR( str(format("Unable to open audit file: %1%") % audit_name) ) ;
  }
  
  void OnStarted( Surface& ) {} 
  
  void OnFinished ( Surface& aSurface )
  { 
    CHECK(aSurface.is_valid());

    CHECK_EQUAL( audit_table.size(), actual_table.size() ) ;
    
    size_t total = audit_table.size() ;
    size_t failed = 0 ;
    for ( size_t i = 0, ei = total ; i != ei ; ++ i )
    {
      size_t idx = i * 2 ;
      Data_ptr  audit_data =  audit_table[idx];
      Data_ptr actual_data = actual_table[idx];

      CHECK(audit_data);
      CHECK(actual_data);

      match cost_m      = equal_cost     (audit_data->cost     ,actual_data->cost);
      match placement_m = equal_placement(audit_data->placement,actual_data->placement);
      
      if ( !cost_m.ok() )
      {
        cerr << "Cost mismatch detected: " << cost_m 
             << "\nExpected: " << audit2str(audit_data) 
             << "\nGot: " << audit2str(actual_data) << endl ;
        ++ failed ;     
      } 
      if ( !placement_m.ok() )
      {
        cerr << "Placement mismatch detected: " << placement_m 
             << "\nExpected: " << audit2str(audit_data) 
             << "\nGot: " << audit2str(actual_data) << endl ;
        ++ failed ;     
      } 
    }
    
    if ( ( failed * 100 / total ) >= 5 )
      throw runtime_error("");
    
  } 
  
  void OnStopConditionReached( Profile const& ) {} 
  
  void OnCollected( Profile const& aProfile, optional<NT> const& aCost, optional<Point> const& aP )
  {
    TRACE( str ( format("I %1% # %2%") % aProfile.v0_v1()->id() % edge2str(aProfile.v0_v1()) ) ) ;
    
    actual_table.insert(make_pair(aProfile.v0_v1()->id(), Data_ptr( new Data(aProfile.v0_v1()->id(),aCost,aP) ) ) ) ;
  }                
  
  void OnSelected( Profile const&, optional<NT> const&, size_t, size_t ) {}
  
  void OnCollapsing( Profile const&, optional<Point> const& ) {}
  
  void OnNonCollapsable( Profile const& ) {}                

  NT toNT ( string s ) 
  { 
    NT r(-1);
    
    try
    {
     r = lexical_cast<NT>(s); 
    }
    catch (...)
    {
      cerr << "Cannot convert string[" << s << "] to a numeric value" << endl ;
    }
    
    return r ;
  }
  
  
  NT compute_epsilon_from_smallest_sample ( NT n ) { return n / NT(256);  }  
  
  struct match
  {
    match ( std::string aFailure ) : mFailure(aFailure) {}
    
    bool ok() const { return mFailure.empty() ; }    
    
    friend std::ostream& operator<< ( std::ostream& os, match const& m )  { return os << m.mFailure ; }
    
    std::string mFailure ;
  } ;
  
  match equal_cost ( optional<NT> const& a, optional<NT> const& b )
  {
    std::string failure ;
    
    if ( a && b )
    { 
      NT actual_diff = CGAL_NTS abs(*a-*b) ;
      
      if ( actual_diff > epsilon_cost )
        failure = str(format("actual_diff=%1% max_diff=%1%") % actual_diff % epsilon_cost) ;
    }
    else if ( !a && b )
    {
      if ( *b < infinite_cost )
      {
        failure = "non-collapsable in case A but collapsable in case B" ;
      }
    }
    else  if ( a && !b )
    {
      if ( *a < infinite_cost )
      {
        failure = "non-collapsable in case B but collapsable in case A" ;
      }
    }
      
    return match(failure);   
  }
  
  match equal_placement ( optional<Point> const& a, optional<Point> const& b )
  {
    std::string failure ;
    
    if ( a && b )
    {
      NT actual_diff = squared_distance(*a,*b) ;
      
      if ( actual_diff > epsilon_sqdist )
        failure = str(format("actual_diff=%1% max_diff=%1%") % actual_diff % epsilon_sqdist) ;
    }
      
    return match(failure);   
  }

  void error ( char const* file, int line, char const* pred, string msg )
  {
    cerr << "ERROR in " << file << " at " << line << endl ;
    if ( pred )
      cerr << "  Assertion failed: " << pred << endl ;
    cerr << "  " << msg << endl ;
   
     
    throw runtime_error("");
  }

private :

  Table audit_table ;
  Table actual_table ;      
  NT    epsilon_cost ;
  NT    epsilon_sqdist ;
  NT    infinite_cost ;
} ;
