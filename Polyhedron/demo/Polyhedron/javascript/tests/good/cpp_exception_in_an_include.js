try {
    include("../bad/uncaught_cpp_exception.js")
} catch(e) {
    print("Exception caught in cpp_exception_in_an_include.js: " + e)
    print("Backtrace:")
    print_backtrace(e.backtrace)
}
