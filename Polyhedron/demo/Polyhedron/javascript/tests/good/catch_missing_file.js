try {
    include("missing-file")
} catch (e)
{
    print_exception_and_bt(e)
    if(!e.toString().match(/No such file or directory/))
        throw "Wrong exception!"
}
print("OK in catch_missing_file")
