try {
    main_window.throw_exception()
} catch(e) {
    print("Exception caught in cpp_exception.js: " + e)
    print("Backtrace:")
    print_backtrace(e.backtrace)
}

print("OK at the end of cpp_exception.js")
