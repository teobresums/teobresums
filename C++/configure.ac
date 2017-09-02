AC_INIT([TEOBResumS], [1.0], [walter.delpozzo@unipi.it])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_LANG_CPLUSPLUS
AC_CONFIG_SRCDIR([src])
AM_INIT_AUTOMAKE
mkdir -p m4
AC_CONFIG_MACRO_DIR([m4])

AM_PROG_AR
AC_PROG_CXX
AC_PROG_LIBTOOL
AC_PROG_INSTALL
LT_INIT
LT_LANG([C++])

AM_CXXFLAGS="-std=c++11 -Ofast -Wall -fPIC"

AC_ARG_ENABLE(debug,
    AC_HELP_STRING([--enable-debug],
    [Add -g to CXXFLAGS]),
    [ CXXFLAGS="-g -O0 $CXXFLAGS"], [ CXXFLAGS="$CXXFLAGS"])

AC_ARG_ENABLE(profile,
              [  --enable-profile        Add gprof-friendly stuff to CXXFLAGS],
              [ CXXFLAGS="$CXXFLAGS -pg -fprofile-arcs" ],)

dnl GSL
AC_LANG_C
AC_ARG_WITH(gsllib, 
            [  --with-gsllib=DIR       DIR=-LDIR, has libgsl and libgslcblas],
            gsllib="$withval/",gsllib= )
if test $gsllib; then
  AC_MSG_RESULT( Checking User specified GSL )
  LDFLAGS="$LDFLAGS -L$gsllib"
  AC_CHECK_LIB( gsl, main, ,
      AC_MSG_ERROR("libgsl not in location specified by --with-gsllib: $gsllib"),
      [ -lm -lgslcblas -L$gsllib ] )
  AC_CHECK_LIB( gslcblas, main, ,
      AC_MSG_ERROR("libgslcblas not in location specified by --with-gsllib: $gsllib"),
      [ -lm -L$gsllib ] )
else
  AC_MSG_RESULT( Checking GSL )
  AC_CHECK_LIB( gsl, main, ,
      AC_MSG_ERROR("libgsl not found; use --with-gsllib"), [-lm -lgslcblas] )
  AC_CHECK_LIB( gslcblas, main, ,
      AC_MSG_ERROR("libgslcblas not found; use --with-gsllib"), [-lm] )
fi

AC_ARG_WITH(gslinc,
        [  --with-gslinc=VAL       VAL=-IVAL, has gsl include files ],
	[gslinc="$withval/" CPPFLAGS="$CPPFLAGS -I$gslinc"],
        [gslinc=] )
AC_MSG_RESULT( Checking GSL Headers )
AC_CHECK_HEADERS([gsl/gsl_errno.h gsl/gsl_math.h gsl/gsl_min.h ],,
  [ if test $gslinc; then
   AC_MSG_ERROR("gsl headers not in location specified by --with-gslinc")
   elif test $gsllib; then
    AC_CHECK_HEADER($gsllib/../include/gsl/gsl_math.h,
     [gslinc="$gsllib/../include/" CPPFLAGS="$CPPFLAGS -I$gslinc"],
      AC_MSG_ERROR("Can't find gsl includes in user specified location; use --with-gslinc"))
    else
      AC_MSG_ERROR("Can't find gsl includes; use --with-gslinc")
   fi ] )


AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([
Makefile
src/Makefile
])

# link any executables against TEOBResumS library
AC_SUBST([LDADD], ["libTEOBResumS.la"])
#AM_CPPFLAGS="-Isrc `gsl-config --cflags` ${AM_CPPFLAGS}"
AC_OUTPUT
AC_MSG_NOTICE([
==================================================
TEOBResumS has now been successfully configured and will be installed under the directory:

${prefix}

Now run 'make' to build TEOBResumS,
and run 'make install' to install TEOBResumS.

==================================================
])
