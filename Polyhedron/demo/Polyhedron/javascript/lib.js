function include(filename) {
    var basename = current_filename.replace(RegExp('/[^/]*$'), "")
    main_window.loadScript(basename + "/" + filename)
}

print_backtrace = function(bt, prefix) {
    var p = typeof prefix !== 'undefined' ?  prefix : "  ";
    for(var i = 0; i < bt.length; ++i) {
        print(p + bt[i])
    }
}

print_exception_and_bt = function(e) {
    print("Caught exception in " + current_filename + ": " + e)
    print("Backtrace:")
    print_backtrace(e.backtrace)
}

quit = main_window.quit

noop = function () {}

noop()
