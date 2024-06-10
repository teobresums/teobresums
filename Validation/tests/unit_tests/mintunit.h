/*
 * mintunit.h
 *	Mint Unit -- a minimal unit testing framework for C
 *
 * Mint Unit has a TAP like output syntax but unlike TAP implementations it
 * stops on the first failing test following the fail-fast fail-hard idea of
 * being forced to fix tests rather than observe broken ones. This code was
 * influenced by http://www.jera.com/techinfo/jtns/jtn002.html
 *
 * Sample usage of mint unit can be seen below:
 *
 *	#include <stdio.h>
 *	#include "mintunit.h"
 *
 *	static char *
 *	test_foo(void)
 *	{
 *		 mt_assert("Test foo function", foo(5) == 5);
 *		 mt_is_null("Test bar function returns NULL", bar(5));
 *		 return NULL;
 *	}
 *
 *	static char *
 *	run_tests(void)
 *	{
 *		mt_init(1);
 *		mt_run_test("Testing foo", test_foo);
 *		return NULL;
 *	}
 *
 *	int
 *	main(int argc, char **argv)
 *	{
 *		char *test_result = run_tests();
 *		if (!test_result)
 *			printf("\n# All tests passed\n");
 *		printf("# Test Suites run: %d\n# Unit Tests run : %d\n", tests_run, units_run);
 *
 *		return test_result == NULL ? 0 : 1;
 *	}
 *
 * Portions Copyright (c) 2016-2017 Daniel Gustafsson <daniel@yesql.se>
 * Portions Copyright (c) 2002 Jera Design
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef _MINTUNIT_H
#define _MINTUNIT_H

#include <stdio.h>
#include <sys/time.h>

unsigned int tests_run = 0;
unsigned int units_run = 0;

/* mt_init - Initalize with the number of suites */
#define mt_init(tests) do { tests_run = units_run = 0; printf("1..%d\n", tests); } while (0)
/* mt_assert - General purpose assertion */
#define mt_assert(message, test) do { units_run++; if (!(test)) return message; } while (0)
/* mt_is_null - Specifically test for NULL */
#define mt_is_null(message, test) do { units_run++; if ((test) != NULL) return message; } while (0)
/* mt_run_test - Run a test suite */
#define mt_run_test(message, test) do { struct timeval stop, start; float msec = 0.0; char *result = NULL; \
	gettimeofday(&start, NULL); result = test(); tests_run++; \
	gettimeofday(&stop, NULL); msec = ((stop.tv_sec - start.tv_sec) * 1000.0f + (stop.tv_usec - start.tv_usec) / 1000.0f); \
	printf("%sok %d - %s (%.3f ms) %s%s\n", (result ? "not " : ""), tests_run, message, msec, result ? "# " : "", result ?: ""); \
	if (result) return result; } while (0)

#endif /* _MINTUNIT_H */
