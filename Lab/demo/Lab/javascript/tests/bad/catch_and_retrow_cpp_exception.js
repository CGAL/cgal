function throw_exception() {
    main_window.throw_exception();
}

f = function() {
    throw_exception()
}

try {
    f()
} catch(e) {
    print_exception_and_bt(e)
    print("Rethrow the exception...")
    throw e
}

print("OK from catch_and_retrow_cpp_exception")
