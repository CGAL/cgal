function throw_exception() {
    main_window.throw_exception();
}

f = function() {
    throw_exception()
}

try {
    f()
} catch(e) {
    print("Caught exception in catch_and_retrow_cpp_exception.js:\n  " + e)
    print("Backtrace:")
    print_backtrace(e.backtrace)
    print("Rethrow the exception...")
    throw e
}

print("OK from catch_and_retrow_cpp_exception")
