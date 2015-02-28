# Even though it's not used at the end, the 1.4.x should be updated before each release.
AC_INIT([BitSeq], [1.10.x])
AC_PREREQ([2.67])

# Determine version.
BS_VERSION=`grep "^Version" DESCRIPTION | sed -e 's/Version:[[[:space:]]]*\([[0-9.]]*\).*/\1/'`
# Set package version define with current version.
AC_DEFINE_UNQUOTED(BS_VERSION, "${BS_VERSION}")
# Substitue DEF_VERSION with version, this is for Makevars.win.
AC_SUBST(DEF_VERSION, [-DBS_VERSION=\\\"${BS_VERSION}\\\"])
VOID_DUMP=[" \" "]
AC_MSG_NOTICE([Version def: $DEF_VERSION])

# Find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
    echo "could not determine R_HOME"
    exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
AC_LANG(C++)
AC_PROG_CXX

# OpenMP magic.
AC_OPENMP
# Substitute OPENMP_CXXFLAGS.
AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS")

# Add SUPPORT_OPENMP pre-processor flag if relevant.
AS_IF([test "${OPENMP_CXXFLAGS}" == "" ],
      AC_MSG_NOTICE([No Open MP support.]),
      AC_DEFINE([SUPPORT_OPENMP]))

# Own defines.
AC_DEFINE([BIOC_BUILD])

# Write config files.
AC_CONFIG_FILES([src/Makevars src/Makevars.win])
AC_OUTPUT
