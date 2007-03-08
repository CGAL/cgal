#include <map>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

class Visitor
{
private :
  
  struct Data
  {
    Data ( size_t i )
      : 
       id(i)
      ,selected(false)
      ,is_collapsable(false)
      ,order( size_t(-1) )
    {}

    static bool match ( Data const& x, Data const& y )
    {  
      bool r = false ;
      if ( x.selected == y.selected )
      {
        if ( x.selected )
        {
          if (    (  x.is_collapsable ==  y.is_collapsable )
               && ( !x.cost           == !y.cost           )
               && ( !x.placement      == !y.placement      )
             )
            r = true ;
        } 
        else r = true ;
      }
      return r ;  
    }

    size_t           id ;
    bool             selected ; 
    bool             is_collapsable ;
    size_t           order ;
    optional<double> cost ;
    optional<Point>  placement ;
  } ;
  typedef shared_ptr<Data> Data_ptr ;
  
  typedef map<int,Data_ptr> Table ;
 
  
private :
  
  typedef char_separator<char> Separator ;
  typedef tokenizer<Separator> Tokenizer ;
  
  // The audit files are generated in a Windows machine, so if this runs on Linux there is a trailing CR that would mess up the tokenizer.
  string normalize_EOL ( string line )
  {
    string::size_type l = line.length();
    string::size_type d = ( l > 0 && line[l-1] == '\r' ) ? 1 : 0 ; 
    return line.substr(0, l-d ) ;
  }

  void ReadAudit( istream& in )
  {
    size_t order = 0 ;
    string line ;
    while ( getline(in,line) )
    {
      if ( line.length() > 1 )
      {
        line = normalize_EOL(line);

        string h = str ( format("AUDIT: %1%") % line ) ;
        history.push_back(h);
        #ifdef TRACE
        std::cerr << h << std::endl ;
        #endif
 
        Tokenizer tk(line,Separator(" "));
        vector<string> tokens(tk.begin(),tk.end());

        CHECK_MSG( tokens.size() > 1, str(format("Invalid audit line, not enough tokens: %1%") % line) ) ;
        char c = tokens[0][0] ;
        switch ( c )
        {
          case 'I' : 
            {
              CHECK_MSG( tokens.size() > 1, str(format("Invalid audit line of type I, id field missing: %1%") % line) ) ;
              size_t id = lexical_cast<size_t>(tokens[1]);

              CHECK_MSG(audit_table.find(id)==audit_table.end()
                       ,str(format("Invalid audit line of type I, id field with duplicate value: %1%") % line)
                       );

              audit_table.insert(make_pair(id, Data_ptr( new Data(id) ) ) ) ;
            }  
            break ;
            
          case 'S' : 
            {
              CHECK_MSG( tokens.size() > 1, str(format("Invalid audit line of type S, id field missing: %1%") % line) ) ;
              size_t id = lexical_cast<size_t>(tokens[1]);
              Data_ptr data = audit_table[id];

              CHECK_MSG(data
                       ,str(format("Invalid audit line of type S, incorrect id field (doesn't match any previous I line): %1%") % line)
                       );
              
              optional<double> cost ;
              if ( tokens.size() > 2 )
                cost = lexical_cast<double>(tokens[2]);

              data->selected = true ;  
              data->cost     = cost ;
              data->order    = order++;
            }  
            break ;
            
          case 'C' : 
            {
              CHECK_MSG( tokens.size() > 1, str(format("Invalid audit line of type C, id field missing: %1%") % line) ) ;
              size_t id = lexical_cast<size_t>(tokens[1]);
              Data_ptr data = audit_table[id];
              
              CHECK_MSG(data
                       ,str(format("Invalid audit line of type C, incorrect id field (doesn't match any previous I line): %1%") % line)
                       );

              optional<Point> placement ;
              if ( tokens.size() > 4 )
              {
                double x = lexical_cast<double>(tokens[2]);
                double y = lexical_cast<double>(tokens[3]);
                double z = lexical_cast<double>(tokens[4]);
                placement = Point(x,y,z);
              }
                
              data->is_collapsable = true ;
              data->placement = placement ;
            }  
            break ;
            
          case 'N' : 
            {
              CHECK_MSG( tokens.size() > 1, str(format("Invalid audit line of type N, id field missing: %1%") % line) ) ;
              size_t id = lexical_cast<size_t>(tokens[1]);
              Data_ptr data = audit_table[id];

              CHECK_MSG(data
                       ,str(format("Invalid audit line of type N, incorrect id field (doesn't match any previous I line): %1%") % line)
                       );
              
              data->is_collapsable = false ;
            }  
            break ;
            
         default :
           REPORT_ERROR( str(format("Invalid audit line: %1%") % line) ) ;
        }
      }  
    }
  }
  
public :
  
  Visitor ( string audit_name ) 
  {
    string h = str ( format("AUDIT FILE: %1%") % audit_name ) ;
    history.push_back(h);
    #ifdef TRACE
    std::cerr << h << std::endl ;
    #endif
    ifstream in(audit_name.c_str());  
    if ( in )
         ReadAudit(in);
    else REPORT_ERROR( str(format("Unable to open audit file: %1%") % audit_name) ) ;
  }
  
  void OnStarted( Surface& ) { order = 0 ; } 
  
  void OnFinished ( Surface& )
  { 
    CHECK_EQUAL( audit_table.size(), actual_table.size() ) ;
    
    size_t total = audit_table.size() ;
    size_t wrong = 0 ;
    for ( size_t i = 0, ei = total ; i != ei ; ++ i )
    {
      size_t idx = i * 2 ;
      Data_ptr  audit_data =  audit_table[idx];
      Data_ptr actual_data = actual_table[idx];

      CHECK(audit_data);
      CHECK(actual_data);

      if ( !Data::match(*audit_data,*actual_data) )
      {
        cerr << "Mismatch detected.\n Expected: " << audit2str(audit_data) << "\n Got: " << audit2str(actual_data) << endl ;
        ++ wrong ;
      } 
    }

    if ( wrong > 0 )
    {
      cerr << wrong << " mismatches out of " << total << " collapses" << endl ;
      if ( wrong > 20 )
      {
        cerr << "TEST FAILED. Too many mismatched collapses." << endl ;
        throw runtime_error("");
      }
    }
  } 
  
  void OnStopConditionReached( Surface& )
  {
  } 
  
  void OnCollected( Halfedge_handle const& aEdge, Surface& )
  {
    string h = str ( format("I %1% # %2%") % aEdge->id() % edge2str(aEdge) ) ;
    history.push_back(h);
    #ifdef TRACE
    std::cerr << h << std::endl ;
    #endif
    actual_table.insert(make_pair(aEdge->id(), Data_ptr( new Data(aEdge->id()) ) ) ) ;
  }                
  
  void OnSelected( Halfedge_handle const& aEdge, Surface&, optional<double> const& aCost, size_t, size_t )
  {
    string h = str ( format("S %1% %2%") % aEdge->id() % opt2str(aCost) ) ;
    history.push_back(h);
    #ifdef TRACE
    std::cerr << h << std::endl ;
    #endif
    actual_table[aEdge->id()]->selected = true ;
    actual_table[aEdge->id()]->cost     = aCost ; 
    actual_table[aEdge->id()]->order    = order ++ ; 
  }                
  
  void OnCollapsing(Halfedge_handle const& aEdge, Surface&, optional<Point> const& aPlacement ) 
  {
    string h = str ( format("C %1%") % aEdge->id() ) ;
    history.push_back(h);
    #ifdef TRACE
    std::cerr << h << std::endl ;
    #endif
    actual_table[aEdge->id()]->placement      = aPlacement ; 
    actual_table[aEdge->id()]->is_collapsable = true ; 
  }                
  
  void OnNonCollapsable(Halfedge_handle const& aEdge, Surface& ) 
  {
    string h = str ( format("N %1%") % aEdge->id() ) ;
    history.push_back(h);
    #ifdef TRACE
    std::cerr << h << std::endl ;
    #endif
    actual_table[aEdge->id()]->is_collapsable = false ; 
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

  Table          audit_table ;
  Table          actual_table ;      
  size_t         order ;
  vector<string> history ;
} ;
