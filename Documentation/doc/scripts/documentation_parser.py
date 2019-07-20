from pyquery import PyQuery as pq
from collections import defaultdict
from sys import argv 
import os.path as op

# if _in is part of args, return true.
def check_type(_in, args):
  if _in in args:
    return True
  else:
    return False

root_path=argv[1]
d = pq(filename=op.join(op.sep, root_path,'index.xml'), parser="xml")
compounds=[p.text() for p in d('compound').items()]
types=[p.attr('kind') for p in d('compound').items()]
type_map = defaultdict(list) #map <type, name>
dict_map = defaultdict(dict)#map <name, map<member type, member name>>
#FOREACH compounds : fill maps
for i in xrange(0,len(compounds)):
  if check_type(types[i], "typedef"):
    types[i]="type"
  name=d('compound').children("name").eq(i).text()
  members=[p.text() for p in d('compound').eq(i).children("member").items()]
  m_types=[p.attr('kind') for p in d('compound').eq(i).children("member").items()]
  if (not check_type(types[i], ['example', 'file', 'dir', 'page', 'group']) and
      not (types[i] == "namespace" and len(members) == 0) and
      not (types[i] == "enum" and len(members) == 0) ):
    if (types[i] == "class"):#check if the class is a concept class
      compound=d('compound').children("name").eq(i).text().replace('_', '__').replace('::', '_1_1')
      filepath='class'+compound+'.xml'
      total_path=op.join(op.sep, root_path,filepath)
      if(op.isfile(total_path)):
        e = pq(filename=total_path, parser="xml")
        compoundnames=[p.text() for p in e('includes').items()]
        
        if(len(compoundnames) > 1 and compoundnames[0].find("Concept") != -1):
          types[i] = 'Concept '+types[i].lower()
    type_map[types[i]].append(name)

    mtype_map = defaultdict(list)# map<member type, member name>

    #FOREACH member :
    for j in xrange(0,len(members)):
      if(check_type(types[i], ['class', 'Concept class'])
         and m_types[j] == "function"):
        m_types[j]="method"
      if(m_types[j] == "typedef"):
        m_types[j]="type"
      mtype_map[m_types[j]].append(members[j])
	  #end FOREACH member
    dict_map[name]=mtype_map
#end FOREACH compound
indent=""


    

#print
#FOREACH type
for btype in type_map:
  out=btype
  if btype.endswith('s'):
    out+='e'
  print out.title()+'s'
  indent+="  "
  #FOREACH name
  for name in type_map[btype]:
    filepath=""
    if check_type(btype, ['class', 'Concept class']):
      filepath='/class'+name.replace('_', '__').replace('::', '_1_1')+'.xml'
    elif btype == 'namespace':
      filepath='/namespace'+name.replace('_', '__').replace('::', '_1_1')+'.xml'
    templates=[]
    if op.isfile(op.join(op.sep, root_path,filepath)):
      f=pq(filename=op.join(op.sep, root_path,filepath), parser="xml")
      templateparams=f("compounddef").children("templateparamlist").eq(0).children("param").items()
      for param in templateparams:
        template_type=""
        template_name=""
        template_defval=""
        template_type=param.children("type").text()
        template_name=param.children("declname").text()
        template_defval=param.children("defval").text()
        complete_template=""
        if not template_type is None:
          complete_template+=template_type+' '
        if not template_name is None:
          complete_template+=template_name
        if not template_defval is None:
          complete_template+=' = '+template_defval
        templates.append(complete_template)        
        if templates==[]:#if no child was found, just take param.text()
         templates=[t.text() for t in param.items()]
    suffix="<"
    #as template got type, defname and declname, name is twice in template. keep only one of them.
    to_remove=[""]
    for template in templates:
      suffix+=template+', '
    if suffix == "<":
      suffix=""
    if suffix.endswith(', '):
      suffix = suffix[:-2]+'>'
    print indent+name+suffix

    indent+="  "
    #FOREACH mtype
    for mtype in (dict_map[name]):
      out=mtype
      if mtype.endswith('s'):
        out+='e'
      print indent+out.title()+'s'
      indent+="  "
      #FOREACH member
      overload_map = defaultdict(int) #contains the number of times a member has appeared (to manage the overloads)
      templateparams=[]
      for member in dict_map[name][mtype]:
        templates=[]
        args="" # will contain the arguments of the methods and functions
        return_type="" #will contain the return type of a function/method
        
        #look for arguments
        if op.isfile(op.join(op.sep, root_path,filepath)):
          f=pq(filename=op.join(op.sep, root_path,filepath), parser="xml")
          index=0
          memberdefs=[m.text() for m in f("memberdef").items()]
          for i in xrange(0,len(memberdefs)):
            member_names=[member_name.text() for member_name in f('memberdef').eq(i).children("name").items()]
            if f('memberdef').eq(i).children("name").text() == member:
              if (index < overload_map[member]):
                index+=1
              elif (index == overload_map[member]):
                if check_type(mtype, ['function', 'method']):
                  args=[f('memberdef').eq(i).children("argsstring").text()]
                  templateparams=f('memberdef').eq(i).children("templateparamlist").children("param").items()
                if check_type(mtype, ['function', 'method', 'type', 'variable']):
                  return_type=[f('memberdef').eq(i).children("type").text()]
                break;
          #end foreach memberdef
        arguments=""
        for arg in args:
          arguments+=arg
        template_types=[]
        template_names=[]
        for param in templateparams:
          template_type=""
          template_name=""
          template_defval=""
          template_type=param.children("type").text()
          template_name=param.children("declname").text()
          template_defval=param.children("defval").text()
          complete_template=""
          if not template_type is None:
            complete_template+=template_type+' '
          if not template_name is None:
            complete_template+=template_name
          if not template_defval is None:
            complete_template+=' = '+template_defval
          templates.append(complete_template)        
          if templates==[]:#if no child was found, just take param.text()
           templates=[t.text() for t in param.items()]

        prefix="template <"
        for template in templates:
          prefix+=template+', '
        if prefix == "template <":
          prefix=""
        if prefix.endswith(', '):
          prefix = prefix[:-2]+'>\n'+indent+"  "
        for definition in return_type:
          prefix+=definition
        if(prefix != ""):
          prefix+=" "
        print indent+prefix+member+arguments
        overload_map[member]+=1
      #END foreach member
      indent=indent[:-2]
    #END foreach mtype
    indent=indent[:-2]
  #END foreach name
  indent=indent[:-2]
#END foreach type
