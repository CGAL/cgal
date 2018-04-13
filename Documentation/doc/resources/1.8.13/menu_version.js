(function() {
  'use strict';

  var url_re =  /doc\.cgal\.org\/(master|latest|(\d\.\d+))\//;
  var url_local =  /.*\/doc_output\//;
  var all_versions = [
    'master',
    'latest',
    '4.12',
    '4.11',
    '4.10',
    '4.9',
    '4.8',
    '4.7',
    '4.6',
    '4.5',
    '4.4',
    '4.3'
  ];

  function build_select(current_version, current_release) {
    var buf = ['<select>'];

    $.each(all_versions, function(id) {
      var version = all_versions[id];
      buf.push('<option value="' + version + '"');
      if (version == current_version) {
console.log(version);
        buf.push(' selected="selected">');
        if (version[0] == '4') {
          buf.push(current_release);
        } else {
          buf.push(title + ' (' + current_release + ')');
        }
      } else {
        buf.push('>' + version);
      }
      buf.push('</option>');
    });

    buf.push('</select>');
    return buf.join('');
  }

  function patch_url(url, new_version) {
    if(url.includes("doc.cgal.org")){  
      return url.replace(url_re, 'https://doc.cgal.org/' + new_version + '/');
    }
    else{
      return url.replace(url_local, 'https://doc.cgal.org/' + new_version + '/');
    }
  }

  function on_switch() {
    var selected = $(this).children('option:selected').attr('value');
    var url = window.location.href,
        new_url = patch_url(url, selected);
    if (new_url != url) {
      window.location.href = new_url;
    }
  }

  $(document).ready(function() {
      var motherNode=document.getElementById("back-nav");
      var node = document.createElement("LI");
      var spanNode = document.createElement("SPAN");
      var titleNode =document.createTextNode("CGAL Version: "); 
      var textNode = document.createTextNode("4.11");
      spanNode.setAttribute("class", "version_menu");
      spanNode.appendChild(textNode);
      node.appendChild(titleNode);
      node.appendChild(spanNode);
      motherNode.insertBefore(node, motherNode.firstChild);
      var match = url_re.exec(window.location.href);
      if (match) {
        var version = match[1];
        var release = '4.11';
        var select = build_select(version, release);
        spanNode.innerHTML=select;
        $('.version_menu').bind('change', on_switch);
      }
      else {
        match = url_local.exec(window.location.href);
        if (match) {
          var version = '4.11';
          var release = '4.11';
          var select = build_select(version, release);
          spanNode.innerHTML=select;
          $('.version_menu select').bind('change', on_switch);
        }
     }
  });
})(); 
