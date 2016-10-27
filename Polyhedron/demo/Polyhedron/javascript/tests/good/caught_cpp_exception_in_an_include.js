try {
    include("../bad/catch_and_retrow_cpp_exception.js")
} catch(e) {
    print("Caught exception in caught_cpp_exception_in_an_include.js:\n  " + e)
    print("Backtrace:")
    print_backtrace(e.backtrace)
}
