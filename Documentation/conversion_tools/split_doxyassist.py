#!/usr/bin/env python
import os
import xml.etree.ElementTree as xml
import re

################################################################################

# A quick-and-dirty script to split up the old doxyassist.xml and 
# place doxygen configuration files in the correct directories.

################################################################################
            
def write_doxy_config (variable_name, the_list, append, file):
    
    if (the_list == []):
        return
    
    if (append):
        identifier = variable_name.ljust(25) +" +=  "
    else:
        identifier = variable_name.ljust(26) +" =  "
    file.write(identifier)
    
    i = 0
    for l in the_list:
        if (i!=0): temp.write("".ljust(30))
        temp.write(l)
        if (i!=len(the_list)-1):
            temp.write(" \ \n")
        i+=1
    temp.write('\n\n')
            
################################################################################

xml.register_namespace("","http://simply-life.net/doxyassist/doxyassist.xsd")
tree = xml.parse('../doxyassist.xml')

# Find every project
for p in tree.iter(tag='{http://simply-life.net/doxyassist/doxyassist.xsd}project'):  

    # Things to create:
    package_name        = ''
    package_inputs      = []
    package_options     = []
    package_references  = []
    package_image_paths = []    
    package_path        = ''    
    
    package_name = p.find('{http://simply-life.net/doxyassist/doxyassist.xsd}name').text
    if package_name == "CGAL": continue

    # - Clean the name so we can find references to it later...
    package_name = package_name.replace(' ','-')
    package_name = package_name.replace('(','-')
    package_name = package_name.replace(')','-')    
    package_name = package_name.replace(',','-')
    
    print ""
    print "Project: "+package_name
       
    doxygen_node = p.find('{http://simply-life.net/doxyassist/doxyassist.xsd}doxygen')    

                
    # Generate consistant input tags for this project.
    for inputs in p.findall('{http://simply-life.net/doxyassist/doxyassist.xsd}input'):
        if ( inputs.text.find("include") != -1 ):
            package_inputs.append(inputs.text)
            continue
        my_input = inputs.text        
        my_input = my_input.replace('.',"")
        input_split =  filter(None, my_input[1:].split('/'))        
        if (input_split[-1] == 'doc'):
            if (my_input[-1] != '/'): my_input += '/'                        
            my_input = '..'+my_input + input_split[0]+'/'
        else:
            if (my_input[-1] != '/'): my_input += '/'
            my_input = '..'+my_input
        package_path = my_input
        print package_path
        package_inputs.append(my_input)
    
    # Get all the tags for this ndoe.
    for e in doxygen_node:
        if 'name' not in e.attrib: continue
        if (e.attrib["name"] == "IMAGE_PATH"):
            for image_items in e:
                if (image_items.text != package_path+'fig' ):
                    package_image_paths.append(image_items.text)
            continue
        if (e.attrib["name"] == "STRIP_FROM_PATH"):            
            continue
            strip_from_path = e.text                        
            
        if (e.attrib["name"] == "STRIP_FROM_INC_PATH"):
            continue
            strip_from_inc_path = e.text
            
        if (e.attrib["name"] == "EXAMPLE_PATH"):
            continue

        if (e.attrib["name"] == "GENERATE_TAGFILE"):
            continue

        if (e.attrib["name"] == "TAGFILES"):
            continue
            
        package_options.append([e.attrib["name"], e.text])    

    
    # Got dependencies from doxygen logs.
    doxy_location = "../log/Doxyfile-CGAL.CGAL."+package_name
    if (os.path.exists(doxy_location)):
        print "Generating references from "+doxy_location
        my_file = open( doxy_location , 'r')
        all_tags=''
        looking_at_tags = False
        for line in my_file:            
            if "TAGFILES" in line and line[0] != '#':
                looking_at_tags=True
            if (looking_at_tags):
                all_tags += line
            if (looking_at_tags) and not ("\\" in line):
                break
                
        # Messy but works..
        all_tags = all_tags.replace(" ", "")    
        all_tags = all_tags.replace("\n", "")            
        all_tags = all_tags.replace("\t", "")                    
        all_tags = all_tags.replace("TAGFILES", "")    
        all_tags = all_tags.replace("=", "")            
        
        for tag in all_tags.split('\\'):
            d = tag.split('/')[2].split('.')[0]
            package_references.append(d)                    
        my_file.close()
    else:
        print "WARNING:" + doxy_location + " doesn't exist"    
                
    # Write the output files for each package.
    path = '../'+package_path
    if os.path.exists(path):

        doxygen_in = os.path.join(path, 'Doxyfile.in')
        temp = open (doxygen_in, 'w')        

        clean_package_inputs = []        
        for this_input in package_inputs:
            if ('include' in this_input):
                new_input = this_input.split('include')[1]            
                new_input = "${CMAKE_SOURCE_DIR}/"+this_input.split('/')[1]+"/include"+new_input
                clean_package_inputs.append(new_input)
            else:
                new_input = this_input.split('doc')[1]            
                new_input = "${CMAKE_SOURCE_DIR}/"+ this_input.split('/')[1] +"/doc"+new_input
                clean_package_inputs.append(new_input)            

        # This will include the default paramters for the configuration.
        temp.write('@INCLUDE = ${CGAL_DOC_PACKAGE_DEFAULTS}\n\n')

        # Inputs for this package
        write_doxy_config ("INPUT", clean_package_inputs,False, temp)
                
        for options in package_options:            
            write_doxy_config(options[0],[options[1]],False, temp)

        write_doxy_config("IMAGE_PATHS", package_image_paths, True, temp)

        temp.close()

        print 'Writing dependencies file'
        dependencies_file = os.path.join(path, 'dependencies')
        temp = open (dependencies_file, 'w')
        for r in package_references:
            temp.write(r+'\n')
        temp.close()

    else:
        print "WARNING: path doesn't exist: "+ package_path
