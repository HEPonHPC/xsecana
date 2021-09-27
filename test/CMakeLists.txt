
set(TESTS test_compile)
list(APPEND TESTS test_hist)
list(APPEND TESTS test_simple_analysis)
list(APPEND TESTS test_simple_efficiency)
list(APPEND TESTS test_simple_flux)
list(APPEND TESTS test_simple_quad_sum)
list(APPEND TESTS test_simple_signal_estimator)
list(APPEND TESTS test_simple_xsec)
list(APPEND TESTS test_systematic)
list(APPEND TESTS test_template_fit_calculator)

foreach(TEST ${TESTS})
	     add_executable("${TEST}" "${TEST}.cc")
	     target_include_directories("${TEST}" PRIVATE
	     					../include
						${EIGEN3_INCLUDE_DIR}
						${ROOT_INCLUDE_DIRS}
	     )
	     target_link_libraries("${TEST}" PRIVATE
				 ${ROOT_LIBRARIES}
				 ROOT::Minuit2
	     )
	     add_test(NAME "${TEST}" COMMAND "${TEST}")
endforeach(TEST)
