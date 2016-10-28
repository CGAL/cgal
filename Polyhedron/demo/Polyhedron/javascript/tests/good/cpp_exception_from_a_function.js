function throw_cpp_exception() {
    main_window.throw_exception();
}

try {
    throw_cpp_exception()
} catch(e) {
    print("Exception caught in cpp_exception_from_a_function.js: " + e)
    print("Backtrace:")
    print_backtrace(e.backtrace)
    if(!e.toString().match(/Exception thrown in MainWindow::throw_exception/))
        throw "Wrong exception!"
}

print("OK at the end of cpp_exception.js")
