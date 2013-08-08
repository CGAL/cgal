#!/usr/bin/python

import argparse
import glob
import codecs
from lxml import etree
import os
from os import path
from pyquery import PyQuery as pq
import re
import shutil
from sys import argv
from sys import stderr

d = pq(filename="Rewrite.html", parser='html', encoding='UTF-8')
u = d.find('ul')[1]

re_new_old_urls=re.compile('<a.*href=\"(.*)\">.*</a>\s+<code>(.*)</code>')
re_old_url=re.compile('<code>(.*)</code>')

f = codecs.open("res.txt", 'w', encoding='utf-8')

def only_html(text):
  return text.split("/")[-1].split("#")[0]
def only_html_and_prefix(text):
  address_split=text.split("/")
  pkg=address_split[-2]
  html=address_split[-1].split("#")[0]
  return pkg+"/"+html

for el in pq(u).find('li'):
  html_text = pq(el).html()
  res = re_new_old_urls.search(html_text)
  if res:
    f.write(only_html_and_prefix(res.group(2))+' '+only_html_and_prefix(res.group(1))+"\n")
    #f.write(only_html(res.group(2))+' '+only_html(res.group(1))+"\n")
    #print "https://cgal.geometryfactory.com/CGAL/test-old-doc-rewrite-rules/"+res.group(2)
    #print res.group(2).split("/")[1], res.group(1).split("/")[1]
  else:
    res=re_old_url.search(html_text)
    if res:
      stderr.write("Warning: "+res.group(1)+" no link found\n")
      f.write("# "+res.group(1)+"\n")
    else:
      stderr.write("Error: "+html_text+"\n");
  


