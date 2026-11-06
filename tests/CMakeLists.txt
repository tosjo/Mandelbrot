add_executable(mandelbrot_tests test_main.cpp)
target_link_libraries(mandelbrot_tests PRIVATE mandelbrot_core mandelbrot_warnings)
add_test(NAME mandelbrot_tests COMMAND mandelbrot_tests)

# Smoke test the command line tool end to end.
add_test(NAME cli_renders_image
    COMMAND mandelbrot --size 64x48 --iterations 200 --antialias 2 --quiet
            --output ${CMAKE_CURRENT_BINARY_DIR}/smoke.png)
add_test(NAME cli_help COMMAND mandelbrot --help)
add_test(NAME cli_rejects_bad_option COMMAND mandelbrot --bogus)
set_tests_properties(cli_rejects_bad_option PROPERTIES WILL_FAIL TRUE)
