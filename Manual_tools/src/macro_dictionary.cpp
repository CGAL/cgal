/**************************************************************************

  macro_dictionary.cpp
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Dictionary of TeX macro definitions.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$

**************************************************************************/

#include <macro_dictionary.h>

#include <iostream>
#include <lexer.h>
#include <input.h>
#include <string_conversion.h>
#include <error.h>
#include <config.h>
#include <list>
#include <string.h>

// #include <hash_map.h>  // get hash_map from mstring.h

using namespace std;

typedef hash_map < string, Macro_item > Macro_dictionary_scope;

class Active_char_scope {
  bool on[256];
public:
   Active_char_scope() {
    memset(on, 0, 256 * sizeof(bool));
  } bool & operator[] (char c) {
    return on[c];
  }
  bool operator[] (char c) const {
    return on[c];
  }
};

typedef
 std::list < Macro_dictionary_scope > Macro_dictionary;
typedef
 std::list < Active_char_scope > Active_char;

Macro_dictionary macro_dictionary(1);
Active_char active_char(1);

void pushMacroScope()
{
  macro_dictionary.push_front(Macro_dictionary_scope());
  active_char.push_front(active_char.front());
}

void popMacroScope()
{
  if (macro_dictionary.size() > 1) {
    macro_dictionary.pop_front();
    active_char.pop_front();
  } else
    printErrorMessage(MacroStackUnderflowError);
}

bool is_active_char(char c)
{
  return active_char.front()[c];
}

void set_all_active_char(char c, bool tag)
{
  Active_char::iterator i = active_char.begin();
  for (; i != active_char.end(); ++i) {
    (*i)[c] = tag;
  }
}

inline void set_active_char(char c, bool tag)
{
  active_char.front()[c] = tag;
}

void localEraseMacro(const string & macro)
{
  if (macro[0] != '\0' && macro[1] == '\0')
    set_all_active_char(macro[0], false);
  Macro_dictionary::iterator i = macro_dictionary.begin();
  while (i != macro_dictionary.end()) {
    Macro_dictionary_scope::iterator j = (*i).find(macro);
    if (j != (*i).end())
      (*i).erase(j);
    ++i;
  }
}

void
insertMacro(const string & macro,
            const string & filename,
            size_t line, const string & body, size_t n_param)
{
  if (macro_def_switch) {
    cerr << "Macro definition: macro `" << macro << "' = `"
        << body << "'." << endl;
  }
  CC_Assert(macro_dictionary.size() > 0);
  if (macro[0] != '\0' && macro[1] == '\0')
    set_active_char(macro[0], true);
  checkMacroOptEnd(macro);
  macro_dictionary.front()[macro] =
      Macro_item(filename, line, body, n_param);
}

void
insertInternalMacro(const string & macro, ExpandFunction fct,
                    size_t n_param)
{
  if (macro_def_switch) {
    cerr << "Internal def.   : macro `" << macro << "'." << endl;
  }
  CC_Assert(macro_dictionary.size() > 0);
  if (macro[0] != '\0' && macro[1] == '\0')
    set_active_char(macro[0], true);
  checkMacroOptEnd(macro);
  macro_dictionary.front()[macro] =
      Macro_item("<internal macro>", 0, fct, n_param);
}

void
insertInternalMacro(const string & macro,
                    const string & body, size_t n_param)
{
  if (macro_def_switch) {
    cerr << "Internal def.   : macro `" << macro << "' = `"
        << body << "'." << endl;
  }
  CC_Assert(macro_dictionary.size() > 0);
  if (macro[0] != '\0' && macro[1] == '\0')
    set_active_char(macro[0], true);
  checkMacroOptEnd(macro);
  macro_dictionary.front()[macro] =
      Macro_item("<internal macro>", 0, body, n_param);
}


void
insertGlobalMacro(const string & macro,
                  const string & filename,
                  size_t line, const string & body, size_t n_param)
{
  if (macro_def_switch) {
    cerr << "Macro global def: macro `" << macro << "' = `"
        << body << "'." << endl;
  }
  CC_Assert(macro_dictionary.size() > 0);
  localEraseMacro(macro);
  if (macro[0] != '\0' && macro[1] == '\0')
    set_all_active_char(macro[0], true);
  checkMacroOptEnd(macro);
  macro_dictionary.back()[macro] =
      Macro_item(filename, line, body, n_param);
}

void
insertInternalGlobalMacro(const string & macro,
                          ExpandFunction fct, size_t n_param)
{
  if (macro_def_switch) {
    cerr << "Internal gdef.  : macro `" << macro << "'." << endl;
  }
  CC_Assert(macro_dictionary.size() > 0);
  localEraseMacro(macro);
  if (macro[0] != '\0' && macro[1] == '\0')
    set_all_active_char(macro[0], true);
  checkMacroOptEnd(macro);
  macro_dictionary.back()[macro] =
      Macro_item("<internal macro>", 0, fct, n_param);
}

void
insertInternalGlobalMacro(const string & macro,
                          const string & body, size_t n_param)
{
  if (macro_def_switch) {
    cerr << "Internal gdef.  : macro `" << macro << "' = `"
        << body << "'." << endl;
  }
  CC_Assert(macro_dictionary.size() > 0);
  localEraseMacro(macro);
  if (macro[0] != '\0' && macro[1] == '\0')
    set_all_active_char(macro[0], true);
  checkMacroOptEnd(macro);
  macro_dictionary.back()[macro] =
      Macro_item("<internal macro>", 0, body, n_param);
}


bool queryMacro(const string & macro, Macro_dictionary_scope::iterator & j)
{
  Macro_dictionary::iterator i = macro_dictionary.begin();
  while (i != macro_dictionary.end()) {
    j = (*i).find(macro);
    if (j != (*i).end())
      return true;
    ++i;
  }
  return false;
}

const Macro_item & fetchMacro(const string & macro)
{
  static Macro_item undefd_macro("undefined", 0, "");
  Macro_dictionary_scope::iterator j;
  if (queryMacro(macro, j))
    return (*j).second;
  undefd_macro.body = macro;
  return undefd_macro;
}

bool definedMacro(const string & macro)
{
  Macro_dictionary_scope::iterator j;
  return queryMacro(macro, j);
}

bool macroIsTrue(const string & macro)
{
  Macro_dictionary_scope::iterator j;
  if (queryMacro(macro, j)) {
    if ((*j).second.fct)
      return false;
    string body = (*j).second.body;
    crop_string(body);
    return (body == "\\lcTrue" || body == "\\ccTrue");
  }
  return false;
}

void eraseMacro(const string & macro)
{
  if (macro_def_switch) {
    cerr << "Macro erasion:    macro `" << macro << "'." << endl;
  }
  localEraseMacro(macro);
}

string expandFirstMacro(string body, bool expand_only_once )
{
  if (macro_exp_switch)
    cerr << "    #X expansion of: `" << body << "'" << endl;
  while (1) {
    while (!body.empty() && isspace(body[0]))   // remove whitespaces
      body.replace(0, 1, "");
    string::const_iterator i = body.begin();
    if (body.empty() || *i != '\\')
      break;
    ++i;
    if (i == body.end())
      break;
    if (isalpha(*i)) {
      ++i;
      while (i != body.end() && isalpha(*i))
        ++i;
    } else {
      if (*i == '\\') {
        ++i;
        if (i == body.end() || *i != '*')
          --i;
      }
    }
    if (!definedMacro(body.substr(0, i - body.begin())))
      break;
    const Macro_item & item = fetchMacro(body.substr(0, i - body.begin()));
    if (item.n_param > 0 || item.n_opt_at_end > 0)
      break;
    // skip whitespaces behind macro
    string::const_iterator j = i;
    while (j != body.end() && isspace(*j) && *j != SEPARATOR)
      ++j;
    if (item.fct)
      body.replace(0, j - body.begin(),
                   item.fct(body.substr(0, i - body.begin()),
                            &body /* dummy */ , 0, 0));
    else
      body.replace(0, j - body.begin(), item.body);
    if( expand_only_once )
      break;
  }
  if (macro_exp_switch)
    cerr << "    results in: `" << body << "'." << endl;
  return body;
}

string
expandMacro(const string & macro,
            const Macro_item & item,
            string parameters[], size_t n_parameters, size_t n_options)
{

  const size_t cache_size = 9;
  static string expand_cache[cache_size];
  static bool cache_valid[cache_size];

  memset(cache_valid, 0, sizeof(bool) * cache_size);

  bool macro_exp_switch2 = macro_exp_switch
      && !(macro == "\\newcommand") && !(macro == "\\newcommand@mom");
  if (macro_exp_switch2) {
    cerr << '`' << macro << "' Expanded using parameters:" << endl;
    for (size_t i = 0; i < n_parameters + n_options; ++i)
      cerr << "    #" << i + 1 << ": `" << parameters[i] << "'" << endl;
    if (!item.fct)
      cerr << "    with body `" << item.body << "'" << endl;
  }
  if (item.fct) {
    if (macro_exp_switch2) {
      string s = item.fct(macro, parameters, n_parameters, n_options);
      cerr << "    internally expanded to: `" << s << "'." << endl;
      return s;
    }
    return item.fct(macro, parameters, n_parameters, n_options);
  }
  string s = item.body;
  string::size_type i = 0;
  bool expand_only_once = false;
  if (!s.empty()) {
    while (i + 1 < s.size()) {  // Expansion loop. At least two characters
      // remain to be inspected.
      if (s[i] == '#') {
        if (s[i + 1] == '#')
          s.replace(i, 1, "");
        else {
          int hash_len = 2;
          int j = i;
          bool glue_tag = false;
          if (s[i + 1] == 'G' && i + 2 < s.size()) {
            glue_tag = true;
            ++i;
            ++hash_len;
          }
          bool expand_tag = false;
          if (s[i + 1] == 'X' && i + 2 < s.size()) {
            expand_tag = true;
            ++i;
            ++hash_len;
          }
          if (s[i + 1] == 'Y' && i + 2 < s.size()) {
            expand_tag = true;
            expand_only_once = true;
            ++i;
            ++hash_len;
          }
          bool skip_tag = false;
          if (s[i + 1] == 'S' && i + 2 < s.size()) {
            skip_tag = true;
            ++i;
            ++hash_len;
          }
          bool length_tag = false;
          if (s[i + 1] == 'L' && i + 2 < s.size()) {
            length_tag = true;
            ++i;
            ++hash_len;
          }
          bool crop_tag = false;
          if (s[i + 1] == 'C' && i + 2 < s.size()) {
            crop_tag = true;
            ++i;
            ++hash_len;
          }
          if ((s[i + 1] >= '1' && s[i + 1] <= '9') ||
              (s[i + 1] >= 'a' && s[i + 1] <= 'z')) {
            size_t index = size_t(s[i + 1]) - size_t('1');
            if (s[i + 1] >= 'a' && s[i + 1] <= 'z')
              index = size_t(s[i + 1]) - size_t('a') + 9;       //'a'==9
            if (index >= n_parameters + n_options) {
              printErrorMessage(ParamIndexError);
              std::cerr << "macro: " << s << " i=" << i << std::endl << std::endl;
              std::cerr << "index: " << index << " n_parameters: " << n_parameters << std::endl;
              ++i;
            } else {
              string repl;
              if (expand_tag) {
                if (glue_tag) {
                  repl = parameters[index];
                  remove_separator(repl);
                  repl = expandFirstMacro(repl,expand_only_once);
                } else {
                  if (index < cache_size && cache_valid[index])
                    repl = expand_cache[index];
                  else {
                    repl = expandFirstMacro(parameters[index],expand_only_once);
                    expand_cache[index] = repl;
                    cache_valid[index] = true;
                  }
                }
              } else {
                repl = parameters[index];
              }
              if (skip_tag) {
                repl.replace(0, 1, "");
              }
              if (crop_tag) {
                crop_string(repl);
              }
              if (length_tag) {
                repl = int_to_string(repl.size());
              }
              if (expand_tag || skip_tag || crop_tag
                  || length_tag || glue_tag) {
                if (glue_tag) {
                  //cerr << "AAA" << convert_quoted_string_seps(repl) << "BBB" << endl;
                  remove_separator(repl);
                  //cerr << "CCC" << convert_quoted_string_seps(repl) << "BBB" << endl;
                  while (j > 0 && s[j - 1] == SEPARATOR) {
                    --j;
                    ++hash_len;
                  }
                  while (j + hash_len < s.size()
                         && s[j + hash_len] == SEPARATOR) {
                    ++hash_len;
                  }
                }
                s.replace(j, hash_len, repl);
                i = j - 1 + repl.size();
              } else {
                s.replace(j, hash_len, repl + SEPARATOR);
                i = j - 1 + repl.size() + 1;
              }
            }
          } else
            printErrorMessage(ParamIndexError);
        }
      }
      ++i;
    }
  }
  if (macro_exp_switch2)
    cerr << "    expanded to: `" << s << "'." << endl;
  return s;
}

void checkMacroOptEnd(const string & macro)
{
  if (macro.size() > 2) {
    int cnt = 0;
    string::const_iterator i = macro.end();
    --i;
    while (i != macro.begin() && *i == 'o') {
      --i;
      ++cnt;
    }
    while (i != macro.begin() && (*i == 'o' || *i == 'm'))
      --i;
    if (*i != '@' || cnt == 0)
      return;
    string basename(macro.begin(), i);
    Macro_dictionary_scope::iterator j;
    if (queryMacro(basename, j)) {
      if (cnt > int ((*j).second.n_opt_at_end))
        (*j).second.n_opt_at_end = cnt;
    } else {
      cerr << endl << "Error: Define macro " << basename
          << " before defining optional arguments " << macro
          << " in `" << in_string->name()
          << " in line " << in_string->line() << "'.";
      if (stack_trace_switch)
        printErrorMessage(MacroUndefinedError);
    }
  }

}


// EOF //
