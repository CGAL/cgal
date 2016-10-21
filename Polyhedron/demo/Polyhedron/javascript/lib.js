include = function (filename) {
    var basename = __FILE__.replace(RegExp('/[^/]*$'), "")
    main_window.loadScript(basename + "/" + filename)
}
