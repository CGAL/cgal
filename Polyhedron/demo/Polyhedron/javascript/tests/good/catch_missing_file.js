try {
    include("missing-file")
} catch (e)
{
    print_exception_and_bt(e)
}
print("OK in catch_missing_file")
