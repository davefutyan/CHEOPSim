
SHELL := /bin/bash

# define this as very first rule to be sure it is exexuted by calling just "make"
all:	make_all


###############################################################################
#      global and most of the case uses definition 
#      if really needed, they can be updated / added in the makecheops file
###############################################################################

# include fist from include directory of the module, then from the installation area
INCLUDES += -Iinclude -I../include -I${CHEOPS_SW}/include 

CXX_CFLAGS += -std=c++11 -fPIC -Wall -Wno-deprecated -DBOOST_LOG_DYN_LINK -DBOOST_TEST_DYN_LINK

EXT_LIBS = ${EXT_LIB_DIRS} -lboost_log -lboost_thread -lboost_date_time -lboost_log_setup \
          -lboost_program_options -lboost_filesystem -lboost_system -lxerces-c -lcfitsio 

CHEOPS_LIBS = -L${CHEOPS_SW}/lib -llogger -lprogram_params -lutilities -lfits_data_model -lfits_dal

COVERAGE_EXCLUDES = 'boost*' '/usr/*' 'Eigen/*' 'c++/*' 'pqxx' 'xercesc/*' 'xsd/*' '*/install/*'

###############################################################################
#     include the module specific definitions
###############################################################################
include makecheops


###############################################################################
#     Get the version of the programs from the RELEASE_VERSION of the
#     sw-product makecheops file if it exist
###############################################################################
VERSION =`grep RELEASE_VERSION ./makecheops 2>/dev/null | sed s/RELEASE_VERSION// | sed s/r_// | sed s/=// | sed s/" "//g; \
          grep RELEASE_VERSION ../makecheops 2>/dev/null | sed s/RELEASE_VERSION// | sed s/r_// | sed s/=// | sed s/" "//g; \
          grep RELEASE_VERSION ../../makecheops 2>/dev/null | sed s/RELEASE_VERSION// | sed s/r_// | sed s/=// | sed s/" "//g; \
          grep RELEASE_VERSION ../../../makecheops 2>/dev/null | sed s/RELEASE_VERSION// | sed s/r_// | sed s/=// | sed s/" "//g`


###############################################################################
#     rule to build the object files. 
#     it is used for every object file
###############################################################################
obj/%.o : src/%.cxx 
	@mkdir -p obj
	g++ -c -g ${CXX_CFLAGS} ${EXT_INC_DIRS} ${INCLUDES} $< -o $(@D)/$(@F)

###############################################################################
#     rule to build the object files with coverage option. They have the
#     extension .oc
#     it is used for every object file
###############################################################################
obj/%.oc : src/%.cxx 
	@mkdir -p obj
	g++ -c -g ${CXX_CFLAGS} -fprofile-arcs -ftest-coverage ${EXT_INC_DIRS} ${INCLUDES} $< -o $(@D)/$(@F)


###############################################################################
#   rules to build the LIB_TARGET1
###############################################################################

ifeq ($(shell uname),Linux)
LIB_EXT = .so
else
LIB_EXT = .dylib
endif

ifdef LIB_TARGET1
FULL_LIB_TARGET1 = lib/lib${LIB_TARGET1}${LIB_EXT}
endif

${FULL_LIB_TARGET1} : $(addprefix obj/, $(LIB_OBJECT1)) 
	@mkdir -p lib
	g++ ${CXX_CFLAGS} -g -shared -dynamiclib -o $@ $^ ${COMPONENT_LIBS} ${CHEOPS_LIBS} ${EXT_LIBS} 


###############################################################################
#   rules to build the LIB_TARGET2
###############################################################################

ifeq ($(shell uname),Linux)
LIB_EXT = .so
else
LIB_EXT = .dylib
endif

ifdef LIB_TARGET2
FULL_LIB_TARGET2 = lib/lib${LIB_TARGET2}${LIB_EXT}
endif

${FULL_LIB_TARGET2} : $(addprefix obj/, $(LIB_OBJECT2)) 
	@mkdir -p lib
	g++ ${CXX_CFLAGS} -g -shared -dynamiclib -o $@ $^ ${COMPONENT_LIBS} ${CHEOPS_LIBS} ${EXT_LIBS} 


###############################################################################
#   rules to build the PYLIB_TARGET1
###############################################################################

ifdef PYLIB_TARGET1
FULL_PYLIB_TARGET1 = lib/lib${PYLIB_TARGET1}.so
endif

${FULL_PYLIB_TARGET1} : $(addprefix obj/, $(PYLIB_OBJECT1)) 
	@mkdir -p lib
	g++ ${CXX_CFLAGS} -g -shared -dynamiclib -o $@ $^ ${PY_COMPONENT_LIBS} ${CHEOPS_LIBS} ${EXT_LIBS} 
	

###############################################################################
# rules to buile EXEC_TARGET1 and EXEC_TARGET2
###############################################################################

ALL_EXEC_TARGETS = $(EXEC_TARGET1) $(EXEC_TARGET2)


ifdef EXEC_TARGET1
FULL_EXEC_TARGET1 = bin/${EXEC_TARGET1}
endif

ifdef EXEC_TARGET2
FULL_EXEC_TARGET2 = bin/${EXEC_TARGET2}
endif

#   create the program EXEC_TARGET1
#   don't do it if we have no source files. This is the case if just
#   a script is defined as EXEC_TARGET1, which should be installe in 
#   ${CHEOPS_SW}/bin
ifdef EXEC_OBJECT1
src/$(EXEC_TARGET1)_init.cxx : makecheops
	sed "s/target_name/$(EXEC_TARGET1)/;s/target_version/$(VERSION)/;s/svnversion/$(shell svnversion)/" ${CHEOPS_SW}/resources/template_init.cxx > $@

# create the EXEC_TARGET1 program 
${FULL_EXEC_TARGET1} : $(addprefix obj/, $(EXEC_OBJECT1)) obj/${EXEC_TARGET1}_init.o
	@mkdir -p bin
	g++ ${CXX_CFLAGS} -o $@ $^ ${COMPONENT_LIBS} ${CHEOPS_LIBS} ${EXT_LIBS}
	@mkdir -p lib
endif

#   create the program EXEC_TARGET2
#   don't do it if we have no source files. This is the case if just
#   a script is defined as EXEC_TARGET2, which should be installe in 
#   ${CHEOPS_SW}/bin
ifdef EXEC_OBJECT2
src/$(EXEC_TARGET2)_init.cxx : makecheops
	sed "s/target_name/$(EXEC_TARGET2)/;s/target_version/$(VERSION)/;s/svnversion/$(shell svnversion)/" ${CHEOPS_SW}/resources/template_init.cxx > $@

# create the EXEC_TARGET2 program 
${FULL_EXEC_TARGET2} : $(addprefix obj/, $(EXEC_OBJECT2)) obj/${EXEC_TARGET2}_init.o
	@mkdir -p bin
	g++ ${CXX_CFLAGS} -o $@ $^ ${COMPONENT_LIBS} ${CHEOPS_LIBS} ${EXT_LIBS} 
	@mkdir -p lib
endif

###############################################################################
#   rules to build the static libraries for the unit tests with coverage option
#   the object files have the extension .oc
###############################################################################
lib/libcoverage.a :  $(addsuffix c, $(addprefix obj/, $(EXEC_OBJECT1))) $(addsuffix c, $(addprefix obj/, $(EXEC_OBJECT2))) $(addsuffix c, $(addprefix obj/, $(LIB_OBJECT1)))
	@if test "X" != "X$^" ; then \
	mkdir -p lib; \
	ar rvs $@ $^;\
	fi;

###############################################################################
#   rules to build the static libraries for the unit tests without coverage option
#   the object files have the extension .o
###############################################################################
lib/libunit_test.a :  $(addprefix obj/, $(EXEC_OBJECT1)) $(addprefix obj/, $(EXEC_OBJECT2)) $(addprefix obj/, $(LIB_OBJECT1))
	@if test "X" != "X$^" ; then \
	mkdir -p lib; \
	ar rvs $@ $^;\
	fi;


###############################################################################
#   rules to build the unit test programs with and withou coverage
###############################################################################

$(foreach  program,$(CXX_UNIT_TESTS),src/$(program)_init.cxx) : makecheops
	sed "s/target_name/$(subst _init.cxx,,$(subst src/,,$@))/;s/target_version/${VERSION}/;s/svnversion/$(shell svnversion)/" ${CHEOPS_SW}/resources/template_init.cxx > $@


$(CXX_UNIT_TESTS) : % : obj/%.oc obj/%_init.oc
	@mkdir -p bin
	g++ ${CXX_CFLAGS}  -o bin/$@ $^  $(wildcard ../lib/lib*.a) $(wildcard ../../lib/lib*.a) \
	${COMPONENT_LIBS} ${CHEOPS_LIBS} ${EXT_LIBS} -lboost_unit_test_framework -lpthread -lgcov




###############################################################################
#     The PHONY targets which are normally used by the user
###############################################################################

FORCE:

.PHONY: unit_test

make_all_targets: ${TARGET}  ${FULL_LIB_TARGET1} ${FULL_LIB_TARGET2} ${FULL_PYLIB_TARGET1} ${FULL_EXEC_TARGET1} ${FULL_EXEC_TARGET2}

make_all:
	@set -e; \
	for i in $(SUBDIRS) .end_of_list; do \
	   ( if test -f $$i/Makefile; then \
	     cd $$i; $(MAKE) $@; \
	     fi; \
	   ) || exit 1; done ;	
	@$(MAKE) make_all_targets 

clean: 
	@set -e; \
	for i in $(SUBDIRS) .end_of_list; do \
	   ( if test -f $$i/Makefile; then \
	     cd $$i; $(MAKE) $@; \
	     fi; \
	   ) || exit 1; done ;
	@if test -f tests/unit_test/Makefile; then \
		cd tests/unit_test; $(MAKE) $@; \
	fi;	
	@if test "X" != "X${EXEC_OBJECT1}" ; then \
		rm -f ${FULL_EXEC_TARGET1}; \
		echo "rm -f ${FULL_EXEC_TARGET1}"; \
	fi;
	@if test "X" != "X${EXEC_OBJECT2}" ; then \
		rm -f ${FULL_EXEC_TARGET2}; \
		echo "rm -f ${FULL_EXEC_TARGET2}"; \
	fi;
	rm -f ${CLEAN} ${FULL_LIB_TARGET1} ${FULL_LIB_TARGET2} ${FULL_PYLIB_TARGET1} \
	      obj/*.o obj/*.oc src/*_init.cxx $(wildcard ../lib/lib*.a) \
	      $(addprefix bin/, ${CXX_UNIT_TESTS}) \
	      obj/*.gcda obj/*.gcno bin/*_init.py
	rm -f .coverage
	rm -f cxxcoverage.xml pycoverage.xml doc/cppcheck.xml doc/pylint.out
	@if [ -d "./tests/unit_test/coverage" ]; then \
	   rm -rf ./tests/unit_test/coverage; \
	   echo "rm -rf ./tests/unit_test/coverage"; \
	fi
	@if test -L resources/job_order_schema.xsd; then \
	  rm -f resources/job_order_schema.xsd; \
	  echo "rm -f resources/job_order_schema.xsd"; \
	fi;  
	@if test -L conf/program_params_schema.xsd; then \
	  rm -f conf/program_params_schema.xsd; \
	  echo "rm -f conf/program_params_schema.xsd"; \
	fi;  

coverage:
	@set -e; \
	for i in $(SUBDIRS) .end_of_list; do \
	   ( if test -f $$i/Makefile; then \
	     echo "Module $$i"; \
	     cd $$i; $(MAKE) $@; \
	     fi; \
	   ) || exit 1; done ;
		
	@if test -f tests/unit_test/Makefile; then \
	    rm -f lib/libunit_test.a; \
	    $(MAKE) lib/libcoverage.a; \
		cd tests/unit_test; $(MAKE) run_coverage; \
	fi;	

run_coverage:
	@for file in $(wildcard ../../obj/*.gcda); do (rm $$file) || exit 1; done ;
	@for file in $(wildcard ./obj/*.gcda);  do (rm $$file) || exit 1; done ;
	@mkdir -p coverage
	@set -e; \
	for t in $(CXX_UNIT_TESTS); do \
	   ( $(MAKE) $$t ) || exit 1; \
	   ( ./bin/$$t $(wildcard ./resources/*_job_order.xml) $(UNIT_TESTS_PROG_PARAMS) \
	   ) || exit 1; done ;	
	@if test "X" != "X${CXX_UNIT_TESTS}" ; then \
		lcov --capture  --directory ../../obj/ --output-file coverage.info; \
		lcov --directory ../../obj/ --output-file coverage.info --remove coverage.info $(COVERAGE_EXCLUDES); \
		genhtml coverage.info --output-directory coverage/cxxcoverage; \
		rm coverage.info; \
		fi;
	@coverage erase;
	@for t in $(PYTHON_UNIT_TESTS); do \
	    set -x; \
	    if test -f python/$${t}.py; then \
	         ( PYTHONPATH=${CHEOPS_SW}/python coverage3 run \
	           -a --source ${CHEOPS_SW}/python/`echo ${PYTHON_COVERAGE} | sed 's/[.]/\//g'` python/$${t}.py ) || exit 1; \
	         fi; \
	    if test -f bin/$${t}.py; then \
	         ( PYTHONPATH=${CHEOPS_SW}/python coverage3 run \
	           -a --source ${CHEOPS_SW}/python/`echo ${PYTHON_COVERAGE} | sed 's/[.]/\//g'` bin/$${t}.py ) || exit 1; \
	         fi; \
	    set +x; \
	    done;
	@if test -f .coverage ; then \
	   coverage3 report -m; \
	   echo "coverage3 html -d ./coverage/pycoverage"; \
	   coverage3 html -d ./coverage/pycoverage; \
	   fi;

coverage_jenkins:
	@set -e; \
	for i in $(SUBDIRS) .end_of_list; do \
	   ( if test -f $$i/Makefile; then \
	     echo "Module $$i"; \
	     cd $$i; $(MAKE) $@; \
	     fi; \
	   ) || exit 1; done ;
	@if test -f tests/unit_test/Makefile; then \
	     rm -f lib/libunit_test.a; \
	    $(MAKE) lib/libcoverage.a; \
	    $(MAKE) coverage_jenkins_level1; \
	fi;     
	-cppcheck --enable=all --inconclusive --xml --xml-version=2 . --force 2> ./doc/cppcheck.xml
	    
	@# Ignore pylint's exit code because pylint returns a non-zero exit code 
	@# even on warnings.
	-find . -name *.py | xargs pylint > ./doc/pylint.out

coverage_jenkins_level1:
	@set -e; \
	if test -f tests/unit_test/Makefile; then \
	   rm -f lib/libunit_test.a; \
	   $(MAKE) lib/libcoverage.a; \
	   cd tests/unit_test; $(MAKE) run_coverage_jenkins; \
	   cd -;\
	   if [ -d "obj" ]; then \
	      echo "gcovr -r `pwd`/.. -e \".*[tT]est.*\" --xml  `pwd`/obj > cxxcoverage.xml"; \
	      gcovr -r `pwd`/.. -e ".*[tT]est.*" --xml `pwd`/obj > cxxcoverage.xml; \
	   fi; \
	   if test -f tests/unit_test/.coverage ; then \
	      set -x; \
	      if [ -d "python" ]; then \
	         module_dir=`pwd`; \
	         cd ${CHEOPS_SW}/python; \
	         coverage3 erase; \
	         coverage3 combine $${module_dir}/tests/unit_test/.coverage; \
	         coverage3 report -m; \
	         coverage3 xml -o $${module_dir}/tests/unit_test/coverage/pycoverage.xml; \
	         coverage3 erase; \
	         cd -; \
	      fi; \
	      set +x; \
	   fi; \
	fi;

run_coverage_jenkins:
	@for file in $(wildcard ../../obj/*.gcda); do (rm $$file) || exit 1; done ;
	@for file in $(wildcard ./obj/*.gcda);  do (rm $$file) || exit 1; done ;
	@mkdir -p coverage
	@set -e; \
	for t in $(CXX_UNIT_TESTS); do \
	   ( $(MAKE) $$t ) || exit 1; \
	   ( ./bin/$$t $(wildcard ./resources/*_job_order.xml) $(UNIT_TESTS_PROG_PARAMS) \
	     --log_format=XML --log_sink=./coverage/unittest_$$t.xml --log_level=all --report_level=no; \
	   ) || exit 1; \
	done ;
	@coverage erase;
	@set -e; \
	set -x; \
	for t in $(PYTHON_UNIT_TESTS); do \
	   if test -f python/$${t}.py; then \
	      ( PYTHONPATH=${CHEOPS_SW}/python coverage3 run \
	        -a --source ${CHEOPS_SW}/python/`echo ${PYTHON_COVERAGE} | sed 's/[.]/\//g'` python/$${t}.py ) || exit 1; \
	      python_tests="$${python_tests} python/$${t}.py"; \
	   fi; \
	   if test -f bin/$${t}.py; then \
	      ( PYTHONPATH=${CHEOPS_SW}/python coverage3 run \
	        -a --source ${CHEOPS_SW}/python/`echo ${PYTHON_COVERAGE} | sed 's/[.]/\//g'` bin/$${t}.py ) || exit 1; \
	      python_tests="$${python_tests} bin/$${t}.py"; \
	   fi; \
	done; \
	if test "X" != "X$${python_tests}" ; then \
	   ( python3 -m "nose" --with-xunit --xunit-file=./coverage/nosetest.xml $${python_tests} ) || exit 1; \
	fi; \
	set +x;

###############################################################################
#    unit tests, functional tests, performance tests, integration tests
#    and all tests together
###############################################################################

unit_test: FORCE
	@set -e; \
	for i in $(SUBDIRS) .end_of_list; do \
	   ( if test -f $$i/Makefile; then \
	     cd $$i; $(MAKE) $@; \
	     fi; \
	   ) || exit 1; done ;
	@if test -f tests/unit_test/Makefile; then \
	    rm -f lib/libcoverage.a; \
	    $(MAKE) lib/libunit_test.a; \
		cd tests/unit_test; $(MAKE) run_unit_test; \
	fi;

run_unit_test:
	@set -e; \
	for t in $(CXX_UNIT_TESTS); do \
	   ( $(MAKE) $$t ) || exit 1; \
	   ( ./bin/$$t $(wildcard ./resources/*_job_order.xml) $(UNIT_TESTS_PROG_PARAMS)) || exit 1; \
	done ;
	@for t in $(PYTHON_UNIT_TESTS); do \
	   (if test -f python/$$t.py; then \
	      echo "python3 -m \"nose\" python/$$t.py"; \
	      ( python3 -m "nose" python/$$t.py) || exit 1; \
	    fi; \
	      if test -f bin/$$t.py; then \
	      echo "python3 -m \"nose\" bin/$$t.py"; \
	      ( python3 -m "nose" bin/$$t.py) || exit 1; \
	   fi; \
	   ) || exit 1; done ;

###############################################################################

functional_test: FORCE
	@set -e; \
	for i in $(SUBDIRS) .end_of_list; do \
	   ( if test -f $$i/Makefile; then \
	     cd $$i; $(MAKE) $@; \
	     fi; \
	   ) || exit 1; done ;
	@if test -f tests/functional/Makefile; then \
		cd tests/functional; $(MAKE) run_functional_test; \
	fi;

run_functional_test:
	@set -e; \
	for t in $(FUNCTIONAL_TESTS); do \
	   ( ./bin/$$t ) || exit 1; done ;	

###############################################################################

integration_test: FORCE
	@set -e; \
	for i in $(SUBDIRS) .end_of_list; do \
	   ( if test -f $$i/Makefile; then \
	     cd $$i; $(MAKE) $@; \
	     fi; \
	   ) || exit 1; done ;
	@if test -f tests/integration/Makefile; then \
		cd tests/integration; $(MAKE) run_integration_test; \
	fi;

run_integration_test:
	@set -e; \
	for t in $(INTEGRATION_TESTS); do \
	   ( ./bin/$$t ) || exit 1; done ;	

###############################################################################

performance_test: FORCE
	@set -e; \
	for i in $(SUBDIRS) .end_of_list; do \
	   ( if test -f $$i/Makefile; then \
	     cd $$i; $(MAKE) $@; \
	     fi; \
	   ) || exit 1; done ;
	@if test -f tests/performance/Makefile; then \
		cd tests/performance; $(MAKE) run_performance_test; \
	fi;

run_performance_test:
	@set -e; \
	for t in $(PERFORMANCE_TESTS); do \
	   ( ./bin/$$t ) || exit 1; done ;	

###############################################################################

tests: unit_test functional_test integration_test performance_test

###############################################################################
#    build the doxygen documentation
###############################################################################

doc: FORCE
	@set -e; \
	for i in $(SUBDIRS) .end_of_list; do \
	   ( if test -f $$i/Makefile; then \
	     cd $$i; $(MAKE) $@; \
	     fi; \
	   ) || exit 1; done ;	
	
	@if test -f doxygen/DoxyfileCxx; then \
	   doxygen doxygen/DoxyfileCxx ; \
	fi;
	@if test -f doxygen/DoxyfilePython; then \
	   doxygen doxygen/DoxyfilePython ; \
	fi;

make_update :
	@set -e; \
	for i in $(SUBDIRS) .end_of_list; do \
	   ( if test -f $$i/Makefile; then \
	     cd $$i; $(MAKE) $@; \
	     fi; \
	   ) || exit 1; done ;
	cp ${CHEOPS_SW}/resources/Makefile .
	@if test -f tests/unit_test/Makefile; then \
		cp ${CHEOPS_SW}/resources/Makefile tests/unit_test/.; \
	fi;
	@if test -f tests/Makefile; then \
		cp ${CHEOPS_SW}/resources/Makefile tests/.; \
	fi;
	@if test -f tests/functional/Makefile; then \
		cp ${CHEOPS_SW}/resources/Makefile tests/functional/.; \
	fi;
	@if test -f tests/integration/Makefile; then \
		cp ${CHEOPS_SW}/resources/Makefile tests/integration/.; \
	fi;
	@if test -f tests/performance/Makefile; then \
		cp ${CHEOPS_SW}/resources/Makefile tests/performance/.; \
	fi;

install: make_all_targets
	@set -e; \
	for i in $(SUBDIRS) .end_of_list; do \
	   ( if test -f $$i/Makefile; then \
	     cd $$i; $(MAKE) $@; \
	     fi; \
	   ) || exit 1; done ;

	@mkdir -p ${CHEOPS_SW}/bin
	@if test "X" != "X${FULL_EXEC_TARGET1}" ; then \
          cp ${FULL_EXEC_TARGET1} ${FULL_EXEC_TARGET2} ${CHEOPS_SW}/bin; \
          echo "cp ${FULL_EXEC_TARGET1} ${FULL_EXEC_TARGET2} ${CHEOPS_SW}/bin"; \
			if test "$(suffix ${EXEC_TARGET1})" = ".py" ; then \
          		sed "s/target_name/${EXEC_TARGET1}/;s/target_version/${VERSION}/;s/svnversion/$(shell svnversion)/" ${CHEOPS_SW}/resources/template_init.py > bin/$(basename ${EXEC_TARGET1})_init.py; \
				cp bin/$(basename ${EXEC_TARGET1})_init.py ${CHEOPS_SW}/bin; \
				echo "cp bin/$(basename ${EXEC_TARGET1})_init.py ${CHEOPS_SW}/bin"; \
          	fi; \
          	if test "$(suffix ${EXEC_TARGET2})" = ".py" ; then \
          		sed "s/target_name/${EXEC_TARGET2}/;s/target_version/${VERSION}/;s/svnversion/$(shell svnversion)/" ${CHEOPS_SW}/resources/template_init.py > bin/$(basename ${EXEC_TARGET2})_init.py; \
				cp bin/$(basename ${EXEC_TARGET2})_init.py ${CHEOPS_SW}/bin; \
				echo "cp bin/$(basename ${EXEC_TARGET2})_init.py ${CHEOPS_SW}/bin"; \
          	fi; \
        fi;

	@mkdir -p ${CHEOPS_SW}/lib
	@if test "X" != "X${LIB_TARGET1}" ; then \
          cp ${FULL_LIB_TARGET1} ${CHEOPS_SW}/lib; \
          echo "cp ${FULL_LIB_TARGET1} ${CHEOPS_SW}/lib"; \
        fi;
	@if test "X" != "X${LIB_TARGET2}" ; then \
          cp ${FULL_LIB_TARGET2} ${CHEOPS_SW}/lib; \
          echo "cp ${FULL_LIB_TARGET2} ${CHEOPS_SW}/lib"; \
        fi;
	@if test "X" != "X${PYLIB_TARGET1}" ; then \
          cp ${FULL_PYLIB_TARGET1} ${CHEOPS_SW}/lib; \
          echo "cp ${FULL_PYLIB_TARGET1} ${CHEOPS_SW}/lib"; \
        fi;        

	@mkdir -p ${CHEOPS_SW}/include
	@if test "X" != "X${INSTALL_INCL}" ; then \
          cp $(addprefix include/, $(INSTALL_INCL)) ${CHEOPS_SW}/include; \
          echo "cp $(addprefix include/, $(INSTALL_INCL)) ${CHEOPS_SW}/include"; \
        fi;

	@mkdir -p ${CHEOPS_SW}/conf
	@if test "X" != "X${EXEC_TARGET1}" ; then \
          cp $(addsuffix .xml ,$(addprefix conf/, $(ALL_EXEC_TARGETS))) ${CHEOPS_SW}/conf; \
          echo "cp $(addsuffix .xml, $(addprefix conf/, $(ALL_EXEC_TARGETS))) ${CHEOPS_SW}/conf"; \
        fi;
	@if test "X" != "X${INSTALL_CONF}" ; then \
          cp $(addprefix conf/, $(INSTALL_CONF)) ${CHEOPS_SW}/conf; \
          echo "cp $(addprefix conf/, $(INSTALL_CONF)) ${CHEOPS_SW}/conf"; \
        fi;

	@mkdir -p ${CHEOPS_SW}/resources
	@if test "X" != "X$(INSTALL_RESOURCES)" ; then \
          cp $(addprefix resources/, $(INSTALL_RESOURCES)) ${CHEOPS_SW}/resources; \
          echo "cp $(addprefix resources/, $(INSTALL_RESOURCES)) ${CHEOPS_SW}/resources"; \
        fi;

	@mkdir -p ${CHEOPS_SW}/python
	@if test "X" != "X$(INSTALL_PYTHON_PACKAGE)" ; then \
		cp -r $(addprefix python/, $(INSTALL_PYTHON_PACKAGE)) ${CHEOPS_SW}/python; \
		echo "cp -r $(addprefix python/, $(INSTALL_PYTHON_PACKAGE)) ${CHEOPS_SW}/python"; \
	fi; 
        
release: FORCE
	@if test "X" != "X${RELEASE_VERSION}"; then \
		${CHEOPS_SW}/bin/release_cheops_package `svn info . | grep "^URL"` $(RELEASE_VERSION); \
	fi;   	     

