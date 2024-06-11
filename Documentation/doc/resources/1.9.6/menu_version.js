(function() {
  'use strict';

  var url_re =  /(cgal\.geometryfactory\.com\/CGAL\/doc\/|doc\.cgal\.org\/)(master|latest|(\d\.\d+|\d\.\d+\.\d+)(-beta\d)?)\//;
  var url_local =  /.*\/doc_output\//;
  var current_version_local = 'master'
  var all_versions = [
      'master',
      '6.0-beta1',
      '5.6.1',
      'latest',
      '5.5.4',
      '5.4.4',
      '5.3.2',
      '5.2.4',
      '5.1.5',
      '5.0.4',
      '4.14.3',
      '4.13.2',
      '4.12.2',
      '4.11.3',
      '4.10.2',
      '4.9.1',
      '4.8.2',
      '4.7',
      '4.6.3',
      '4.5.2',
      '4.4',
      '4.3'
  ];

  function build_select(current_version) {
    if( current_version == 'master') {
      let top_elt = document.getElementById("top");

      let first_element = top_elt.childNodes[0];
      let new_div = document.createElement("p");
      new_div.innerHTML = '⚠️ This documentation corresponds to the <a style="font-familly: monospace;" href="https://github.com/CGAL/cgal/tree/master">master</a> development branch of CGAL. It might diverge from the official releases.';
      new_div.style.cssText = "background-color: #ff9800; margin: 1ex auto 1ex 1em; padding: 1ex; border-radius: 1ex; display: inline-block;"
      let OK = top_elt.insertBefore(new_div, first_element);
    }
    var buf = ['<select>'];
    $.each(all_versions, function(id) {
      var version = all_versions[id];
      buf.push('<option value="' + version + '"');
      if (version == current_version) {
        buf.push(' selected="selected">' + version);
      } else {
        buf.push('>' + version);
      }
      buf.push('</option>');
    });
    if ( !all_versions.includes(current_version)) {
       buf.push('<option value="' + current_version + '"');
       buf.push(' selected="selected">' + current_version);
       buf.push('</option>');
    }
    buf.push('</select>');
    return buf.join('');
  }

  function patch_url(url, new_version) {
    if(url.includes("doc.cgal.org")||url.includes("cgal.geometryfactory.com")){
      return url.replace(url_re, 'doc.cgal.org/' + new_version + '/');
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
      var motherNode=$("#back-nav ul")[0];
      var node = document.createElement("LI");
      var spanNode = document.createElement("SPAN");
      var titleNode =document.createTextNode("CGAL Version: ");
      var textNode = document.createTextNode("x.y");
      spanNode.setAttribute("class", "version_menu");
      spanNode.appendChild(textNode);
      node.appendChild(titleNode);
      node.appendChild(spanNode);
      motherNode.insertBefore(node, motherNode.firstChild);
      $("#back-nav").css("padding-top", "0").css("padding-bottom", "0");
      var match = url_re.exec(window.location.href);
      if (match) {
        var version = match[2];
        var select = build_select(version);
        spanNode.innerHTML=select;
        $('.version_menu select').bind('change', on_switch);
      }
      else {
        match = url_local.exec(window.location.href);
        if (match) {
          var version = current_version_local;
          var select = build_select(version);
          spanNode.innerHTML=select;
          $('.version_menu select').bind('change', on_switch);
        }
     }
  });
})();
