#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([FULL-PACKAGE-NAME], [VERSION], [BUG-REPORT-ADDRESS])
AM_INIT_AUTOMAKE(hic_breakfinder,1.0)
AC_CONFIG_SRCDIR([src/hic_breakfinder.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lbamtools':
AC_CHECK_LIB([bamtools], [main])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_CHECK_FUNCS([floor pow])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
