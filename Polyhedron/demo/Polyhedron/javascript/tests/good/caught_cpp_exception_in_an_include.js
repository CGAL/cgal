var good = false
try {
    include("../bad/catch_and_retrow_cpp_exception.js")
} catch(e) {
    print("Caught exception in caught_cpp_exception_in_an_include.js:\n  " + e)
    print("Backtrace:")
    print_backtrace(e.backtrace)
    if(!e.toString().match(/Exception thrown in MainWindow::throw_exception/))
        throw "Wrong exception!"
    good = true
}
if(!good) throw "The exception was not caugth!"
print("OK in caught_cpp_exception_in_an_include.js")
