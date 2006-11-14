#include <map>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

class Visitor
{
private :
  
  struct Data
  {
    Data ()
      : 
       selected(false)
      ,is_collapsable(false)
      ,order( size_t(-1) )
    {}

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
  
  void ReadAudit( istream& in )
  {
    size_t order = 0 ;
    string line ;
    while ( getline(in,line) )
    {
      if ( line.length() > 1 )
      {
        //std::cerr << "AUDIT: " << line << " (order=" << order << ")" << std::endl ;
 
        Tokenizer tk(line,Separator(" "));
        vector<string> tokens(tk.begin(),tk.end());
        CHECK( tokens.size() > 1 ) ;
        char c = tokens[0][0] ;
        switch ( c )
        {
          case 'I' : 
            {
              CHECK( tokens.size() > 1 ) ;
              size_t id        = lexical_cast<size_t>(tokens[1]);
              CHECK(audit_table.find(id)==audit_table.end());
              audit_table.insert(make_pair(id, Data_ptr( new Data() ) ) ) ;
            }  
            break ;
            
          case 'S' : 
            {
              size_t id = lexical_cast<size_t>(tokens[1]);
              Data_ptr data = audit_table[id];
              CHECK(data);
              
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
              size_t id = lexical_cast<size_t>(tokens[1]);
              Data_ptr data = audit_table[id];
              CHECK(data);
              
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
              size_t id = lexical_cast<size_t>(tokens[1]);
              Data_ptr data = audit_table[id];
              CHECK(data);
              
              data->is_collapsable = false ;
            }  
            break ;
            
         default :
           SHOW_ERROR("Invalid audit line: " << line );
        }
      }  
    }
  }
  
public :
  
  Visitor ( string audit_name ) 
  {
    //std::cerr << "audit_name: " << audit_name << std::endl ;
    ifstream in(audit_name.c_str());  
    if ( in )
         ReadAudit(in);
    else SHOW_ERROR("Unable to open audit file: " << audit_name);
  }
  
  void OnStarted( Surface& ) { order = 0 ; } 
  
  void OnFinished ( Surface& )
  { 
    CHECK( audit_table.size() == actual_table.size() );
    
    for ( size_t i = 0, ei = audit_table.size() ; i != ei ; ++ i )
    {
      size_t idx = i * 2 ;
      Data_ptr audit_data  = audit_table[idx];
      Data_ptr actual_data = actual_table[idx];
      CHECK(audit_data);
      CHECK(actual_data);
      CHECK_EQUAL(audit_data->selected,actual_data->selected);
      if ( audit_data->selected )
      {
        //CHECK_EQUAL(audit_data->order          ,actual_data->order); 
        CHECK_EQUAL(audit_data->is_collapsable ,actual_data->is_collapsable);
        CHECK_EQUAL(!audit_data->cost          ,!actual_data->cost);
        CHECK_EQUAL(!audit_data->placement     ,!actual_data->placement);
      }
    }
  } 
  
  void OnStopConditionReached( Surface& )
  {
  } 
  
  void OnCollected( Halfedge_handle const& aEdge, Surface& )
  {
    //std::cerr << "ACTUAL: I " << edge2str(aEdge) << std::endl ;
    actual_table.insert(make_pair(aEdge->id(), Data_ptr( new Data() ) ) ) ;
  }                
  
  void OnSelected( Halfedge_handle const& aEdge, Surface&, optional<double> const& aCost, size_t, size_t )
  {
    //std::cerr << "ACTUAL: S " << edge2str(aEdge) << " cost:" << aCost << " (order=" << order << ")" << std::endl ;
    actual_table[aEdge->id()]->selected = true ;
    actual_table[aEdge->id()]->cost     = aCost ; 
    actual_table[aEdge->id()]->order    = order ++ ; 
  }                
  
  void OnCollapsing(Halfedge_handle const& aEdge, Surface&, optional<Point> const& aPlacement ) 
  {
    //std::cerr << "ACTUAL: C"  << aEdge->id() << std::endl ;
    actual_table[aEdge->id()]->placement      = aPlacement ; 
    actual_table[aEdge->id()]->is_collapsable = true ; 
  }                
  
  void OnNonCollapsable(Halfedge_handle const& aEdge, Surface& ) 
  {
    //std::cerr << "ACTUAL: N" << aEdge->id() << std::endl ;
    actual_table[aEdge->id()]->is_collapsable = false ; 
  }                
  
private :

  Table  audit_table ;
  Table  actual_table ;      
  size_t order ;
} ;
