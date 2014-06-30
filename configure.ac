#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(partest, 1.0, ddarriba@udc.es)
AM_INIT_AUTOMAKE(partest, 1.0)

AC_CONFIG_SRCDIR(src/PartitionTest.cpp)
AC_CONFIG_HEADERS(config.h)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LN_S
AC_PROG_CPP

# Checks for libraries.
AC_CHECK_LIB([m], [main])
AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_LIB([pll-sse3], [pllCreateInstance])

#AC_SEARCH_LIBS([pllCreateInstance], [pll-sse3], [have_pllsse=yes])
#if test "x${have_pllsse}" = "xyes"; then
  AC_CHECK_HEADERS([pll.h], [], [have_pllsse=no])
  AM_LDFLAGS=-lpll-sse3
#fi

# Checks for header files.
AC_FUNC_ALLOCA
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor mkdir pow sqrt strrchr strtol])

# AC_CONFIG_FILES([Makefile])
AC_OUTPUT(Makefile src/Makefile man/Makefile)
