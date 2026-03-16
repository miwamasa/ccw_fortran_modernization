# オリジナル VA10A (minimize) テスト
add_executable(test_minimize test_minimize.f minimize.f)
add_test(NAME minimize_test COMMAND test_minimize)

# Levenberg-Marquardt法テスト
add_executable(test_lm lm_solver.f90 test_lm.f90)
add_test(NAME lm_test COMMAND test_lm)

# Broyden法テスト
add_executable(test_broyden broyden_solver.f90 test_broyden.f90)
add_test(NAME broyden_test COMMAND test_broyden)
