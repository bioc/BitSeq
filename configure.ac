m4_define([BS_VERSION],
          m4_bpatsubst(m4_esyscmd([grep Version DESCRIPTION]),
                       [.*Version:[ \t]\([0-9.]*\)\(\w\|\W\)*],
                       [\1]))
AC_INIT([Bitseq], [BS_VERSION], [glaus@cs.man.ac.uk])
AC_PREREQ([2.59])

m4_include([m4/ac_openmp.m4])

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
AC_PROG_CXX

AC_LANG(C++)

# OpenMP magic.
AC_OPENMP
AC_SUBST(OPENMP_CXXFLAGS,${OPENMP_CXXFLAGS})

# Add SUPPORT_OPENMP pre-processor flag if relevant.
AS_IF([test "${OPENMP_CXXFLAGS}" == "" ],
      AC_MSG_NOTICE([No Open MP support.]),
      AC_DEFINE([SUPPORT_OPENMP]))

# Own defines.
AC_DEFINE([BIOC_BUILD])

# Rsamtools makefile includes:
AC_SUBST(SAM_INCLUDE, [`echo 'cat(system.file("usretc", .Platform[["r_arch"]],"Rsamtools.mk", package="Rsamtools", mustWork=TRUE))' |"${R_HOME}/bin/R" --vanilla --slave`])
AC_MSG_NOTICE([SAM_INCLUDE: $SAM_INCLUDE])

# Set package version substitution for windows makefile.
AC_SUBST(PACKAGE_VERSION, ${PACKAGE_VERSION})

# Write config files.
AC_CONFIG_FILES([src/Makevars src/Makevars.win])
AC_OUTPUT
