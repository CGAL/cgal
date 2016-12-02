var good = false
try {
    include("missing-file")
} catch (e)
{
    print_exception_and_bt(e)
    if(!e.toString().match(/No such file or directory/))
        throw "Wrong exception!"
    good = true
}
if(!good) throw "The exception was not caugth!"
print("OK in catch_missing_file")
