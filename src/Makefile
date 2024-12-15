OS_VERSION=$(shell uname -s)

VERTER_FLAG=-c -g
S21_FLAGS=-c -g -Wall -Werror -Wextra -std=c11 
ALL_SOURCES=$(shell find . -name "*.[c|h]")
PRE_TEST_FILE=matrix_test
TEST_FILE=matrix_test.c

COVERAGE=--coverage


s21_matrix.a:
	gcc $(VERTER_FLAG) s21_matrix.c -o s21_matrix.o
	ar rcs s21_matrix.a s21_matrix.o
	ranlib s21_matrix.a


all: test_build cppcheck style gcov_report leaks

test_build:
	gcc -g $(S21_FLAGS) s21_matrix.c $(COVERAGE) -o libs21_matrix.o
	ar rcs libs21_matrix.a libs21_matrix.o
	ranlib libs21_matrix.a

test.c:
	-rm $(TEST_FILE)
	checkmk $(PRE_TEST_FILE) > $(TEST_FILE)

test: test.c test_build
	gcc matrix_test.c -lcheck -L. -ls21_matrix -lgcov 
	./a.out

test_ub: test.c test_build
	gcc matrix_test.c -g -lcheck -lm -lrt -lsubunit -pthread -L . -ls21_matrix -lgcov
	./a.out

gcov_report: test
	lcov -t "s21_matrix" -o s21_matrix.info -c -d .
	genhtml -o report s21_matrix.info
	open report/index.html

gcov_ub: test_ub test_build
	lcov -t "s21_matrix" -o s21_matrix.info -c -d .
	genhtml -o report s21_matrix.info
	open report/index.html

valgrind: test_ub
	valgrind --leak-check=full --show-leak-kinds=all --tool=memcheck --trace-children=yes --track-origins=yes ./a.out

leaks: test
	CK_FORK=no leaks -atExit -- ./a.out

cppcheck: test.c
	cppcheck --enable=all --suppress=missingIncludeSystem *.c *.h

style:
	clang-format -n -style=google $(ALL_SOURCES)

fix_style:
	clang-format -i -style=google $(ALL_SOURCES)

clean:
	-rm *.o *.out *.a *.gcno *.gcda *.info matrix_test.c
	-rm -rf report/
add: fix_style
	git add $(ALL_SOURCES) $(PRE_TEST_FILE) Makefile



