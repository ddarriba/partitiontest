#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.50])
AC_INIT(partest, 1.0, ddarriba@udc.es)
AM_INIT_AUTOMAKE # (partest, 1.0)

AC_CONFIG_SRCDIR([src/PartitionTest.cpp])
AC_CONFIG_HEADERS([config.h])

AC_LANG(C++)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LN_S
AC_PROG_CPP

# Checks for libraries.
AC_CHECK_LIB([m], [main])
AC_CHECK_LIB([pthread], [pthread_create])

AC_SEARCH_LIBS([pllCreateInstance], [pll-avx-pthreads], [have_pthreads=yes have_avx=yes], [have_pthreads=no])
if test "x${have_pthreads}" = "xno"; then
  AC_SEARCH_LIBS([pllCreateInstance], [pll-sse3-pthreads], [have_pthreads=yes have_sse=yes], [have_pthreads=no])
  if test "x${have_pthreads}" = "xno"; then
    AC_SEARCH_LIBS([pllCreateInstance], [pll-avx], [have_avx=yes], [have_avx=no])
    if test "x${have_avx}" = "xno"; then
      AC_SEARCH_LIBS([pllCreateInstance], [pll-sse3], [have_sse=yes], [have_sse=no])
    fi
  fi
fi
AM_CONDITIONAL(HAVE_PTHREADS, test "x$have_pthreads" = "xyes")
AM_CONDITIONAL(HAVE_AVX, test "x$have_avx" = "xyes")
AM_CONDITIONAL(HAVE_SSE, test "x$have_sse" = "xyes")

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

have_mpi=no
AX_MPI([have_mpi=yes],[have_mpi=no
  if test x"$have_mpi" = xyes; then
     AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
  else
     AC_MSG_WARN([No MPI compiler found, won't use MPI.])
  fi
])
AC_SUBST(MPICXX)

AM_CONDITIONAL(BUILD_MPI, test "x$have_mpi" = "xyes")

AC_CONFIG_FILES([Makefile
                 man/Makefile
                 mpi/Makefile
                 src/Makefile])
AC_OUTPUT
