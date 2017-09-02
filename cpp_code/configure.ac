mkdir -p m4
for f in '-std=c++11' '-std=c++11 -stdlib=libc++'
do
AX_CHECK_COMPILE_FLAG([$f], [CXXFLAGS="$CXXFLAGS $f" stdpass=true], [], [], [])
${stdpass-false} && break
done
if ! "${stdpass-false}"; then
AC_MSG_ERROR([Unable to turn on C++11 mode with this compiler])
fi
AM_CXXFLAGS="-std=c++11 -Ofast -Wall -fPIC"

AC_INIT([TEOBResum], [1.0], [walter.delpozzo@unipi.it])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_LANG([C])
AC_LANG([C++])
AC_PROG_CXX
AM_PROG_AR
LT_INIT
LT_LANG([C++])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([
Makefile
src/Makefile
])
AC_CONFIG_MACRO_DIRS([m4])
ACLOCAL_AMFLAGS+="-I m4"
# check for required libraries
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
# link any executables against TEOBResum library
AC_SUBST([LDADD], ["libTEOBResum.la"])
AM_CPPFLAGS="-Isrc `gsl-config --cflags` ${AM_CPPFLAGS}"
AC_OUTPUT
AC_MSG_NOTICE([
==================================================
TEOBResum has now been successfully configured:

run 'make' to build TEOBResum,
and run 'make install' to install TEOBResum.

==================================================
])
