var good = false
try {
    include("../bad/uncaught_cpp_exception.js")
} catch(e) {
    print("Exception caught in cpp_exception_in_an_include.js: " + e)
    print("Backtrace:")
    print_backtrace(e.backtrace)
    if(!e.toString().match(/Exception thrown in MainWindow::throw_exception/))
        throw "Wrong exception!"
    good = true
}
if(!good) throw "The exception was not caugth!"
print("OK in cpp_exception_in_an_include.js")
