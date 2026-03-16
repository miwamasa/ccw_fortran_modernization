add_executable(test_minimize test_minimize.f minimize.f)

add_test(NAME minimize_test COMMAND test_minimize)
