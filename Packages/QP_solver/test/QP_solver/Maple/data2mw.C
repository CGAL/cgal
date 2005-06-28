#include <iostream>
#include <fstream>
#include <string>

const std::string Header = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Worksheet><Version major=\"6\" minor=\"1\"/><View-Properties><Zoom percentage=\"100\"/></View-Properties><Styles><Layout alignment=\"left\" bullet=\"none\" linespacing=\"0.0\" name=\"Heading 1\" spaceabove=\"8.0\" spacebelow=\"4.0\"/><Layout alignment=\"left\" bullet=\"none\" firstindent=\"0.0\" leftmargin=\"0.0\" linebreak=\"space\" linespacing=\"0.0\" name=\"Normal\" rightmargin=\"0.0\" spaceabove=\"0.0\" spacebelow=\"0.0\"/><Font background=\"[0,0,0]\" bold=\"true\" executable=\"true\" family=\"Monospaced\" foreground=\"[255,0,0]\" name=\"Maple Input\" opaque=\"false\" size=\"12\"/><Font background=\"[0,0,0]\" bold=\"false\" executable=\"false\" family=\"Times New Roman\" foreground=\"[0,0,0]\" italic=\"false\" name=\"Text\" opaque=\"false\" size=\"12\" underline=\"false\"/><Font background=\"[0,0,0]\" bold=\"true\" family=\"Serif\" name=\"Heading 1\" opaque=\"false\" size=\"18\"/></Styles>";
const std::string Footer = "<Text-field/></Worksheet>\n";
const std::string  Pattern1 = "%%%%%%%1";

const std::string Section = "<Section><Title><Text-field layout=\"Heading 1\" style=\"Heading 1\">%%%%%%%1</Text-field></Title>";
const std::string SectionEnd = "</Section>";

const std::string InputGroup = "<Group><Input><Text-field layout=\"Normal\" prompt=\"&gt; \" style=\"Maple Input\">";
const std::string TextGroup = "<Group><Input><Text-field layout=\"Normal\" style=\"Text\">";
const std::string GroupEnd = "</Text-field></Input></Group>";

const std::string SectionToken = "### Section: ";
const std::string MapleCommentToken = "## ";

//#define DEBUG(x) { x; }
#define DEBUG(x) ;

void bailout(bool assertion,char *msg)
{
  if (!assertion) {
    std::cerr << msg << std::endl;
    exit(2);
  }
}

std::string unspacify(const std::string& s,bool fromRight=true)
{
  // find starting position:
  unsigned int i = 0;
  while (i < s.length() && s[i] == ' ')
    ++i;
  
  // find ending position:
  int j = s.length()-1;
  if (fromRight)
    while (j > 0 && s[j] == ' ')
      --j;
  return s.substr(i,j-i+1);
}

std::string replace1(std::string s,const std::string arg1)
{
  return s.replace(s.find(Pattern1),Pattern1.length(),arg1);
}

int main(int argn,char **argv)
{
  std::ifstream inf;
  std::ofstream outf;

  // check input:
  bailout(argn == 3,"Error: usage is 'data2mw input-file output-file");
  inf.open(argv[1]);
  bailout(!(!inf),"Error: could not open input file");
  outf.open(argv[2]);
  bailout(!outf.fail(),"Error: could not open output file");
  outf << Header;

  bool isSectionOpen = false;
  bool isInputGroupOpen = false;
  bool isTextGroupOpen = false;
  std::string line, content;
  while(inf.good()) { 
    bailout(!(isInputGroupOpen == true && isTextGroupOpen == true),
	    "Error: internal error");
    std::getline(inf,line);
    line = unspacify(line,false);
    DEBUG(std::cout << "LINE IS: " << line << std::endl;);
    
    // empty lines end groups:
    if (line.length() == 0) {
      if (isTextGroupOpen || isInputGroupOpen)
	outf << GroupEnd;
      isTextGroupOpen = isInputGroupOpen = false;
    }
    
    // nonempty lines are parsed:
    else {
      if (line.find(SectionToken) == 0) {
	content = unspacify(line.substr(SectionToken.length()));
	// close groups:
	if (isTextGroupOpen || isInputGroupOpen)
	  outf << GroupEnd;
	isTextGroupOpen = isInputGroupOpen = false;
	// close section:
	if (isSectionOpen)
	  outf << SectionEnd;
	
	DEBUG(std::cout << "CONTENT IS: " << content << std::endl
                        << "RESULT IS: " << replace1(Section,content)
                        << std::endl;);
	outf << replace1(Section,content);
	isSectionOpen = true;
      } else if (line.find(MapleCommentToken) == 0) {
	content = unspacify(line.substr(MapleCommentToken.length()));
	if (isInputGroupOpen) {
	  outf << GroupEnd;
	  isInputGroupOpen = false;
	}
	if (!isTextGroupOpen) {
	  outf << TextGroup;
	  isTextGroupOpen = true;
	} else
	  outf << std::endl;
	outf << content;
      } else {
	if (isTextGroupOpen) {
	  outf << GroupEnd;
	  isTextGroupOpen = false;
	}
	if (!isInputGroupOpen) {
	  outf << InputGroup;
	  isInputGroupOpen = true;
	} else
	  outf << std::endl;
	outf << line;
      }
    }
  }
  
  // close groups and section:
  if (isTextGroupOpen || isInputGroupOpen)
    outf << GroupEnd;
  if (isSectionOpen)
    outf << SectionEnd;
  outf << Footer;
  inf.close();
  outf.close();
}
