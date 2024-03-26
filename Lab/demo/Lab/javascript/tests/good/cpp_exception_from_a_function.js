var good = false
function throw_cpp_exception() {
    main_window.throw_exception();
}

try {
    throw_cpp_exception()
} catch(e) {
    print_exception_and_bt(e)
    if(!e.toString().match(/Exception thrown in MainWindow::throw_exception/))
        throw "Wrong exception!"
    good = true
}
if(!good) throw "The exception was not caught!"

print("OK at the end of cpp_exception.js")
