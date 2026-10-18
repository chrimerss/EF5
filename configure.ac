AC_INIT([Ensemble Framework For Flash Flood Forecasting], [0.1], [zac.flamig@noaa.gov],
             [ef5], [http://hydro.ou.edu/])
AC_PREREQ([2.62])
AC_PROG_CXX
AC_OPENMP
AC_SUBST(OPENMP_CFLAGS)

# C++11 is required (threads in GridWriterFull). Compilers before GCC 6
# default to C++98, so ask for C++11 when the default is older.
AC_LANG_PUSH([C++])
m4_define([EF5_CXX11_TEST], [AC_LANG_PROGRAM([[#if __cplusplus < 201103L
#error C++11 required
#endif]], [])])
AC_MSG_CHECKING([whether $CXX supports C++11 by default])
AC_COMPILE_IFELSE([EF5_CXX11_TEST],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   CXX="$CXX -std=gnu++11"
   AC_MSG_CHECKING([whether $CXX works])
   AC_COMPILE_IFELSE([EF5_CXX11_TEST],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      AC_MSG_ERROR([EF5 needs a C++11 compiler])])])
AC_LANG_POP([C++])

AC_CANONICAL_HOST
case $host_os in
    *mingw*)
        CFLAGS="-static -static-libgcc -I/usr/x86_64-w64-mingw32/sys-root/mingw/include $CFLAGS"
        CXXFLAGS="-static -static-libgcc -I/usr/x86_64-w64-mingw32/sys-root/mingw/include $CXXFLAGS"
        CPPFLAGS="-static -static-libgcc -I/usr/x86_64-w64-mingw32/sys-root/mingw/include $CPPFLAGS"
        LDFLAGS="-static -static-libgcc -L/usr/x86_64-w64-mingw32/sys-root/mingw/lib $LDFLAGS"
	WINDOWS=yes
	AC_CHECK_TOOL(WINDRES, windres, no)	
        ;;
     *)
	;;
esac
AM_CONDITIONAL(WINDOWS, test x$WINDOWS = xyes)
#AC_CHECK_LIB([z], [gzread], [],[
#         echo "z library is required for this EF5"
#         exit -1])
#AC_CHECK_LIB([tiff], [TIFFReadScanline], [],[
#         echo "Tiff library is required for this EF5"
#         exit -1])
#AC_CHECK_LIB([geotiff], [GTIFNew], [],[
#         echo "GeoTiff library is required for this EF5"
#         exit -1])
#if test x$WINDOWS != xyes; then
#AC_CHECK_LIB([gomp], [omp_get_wtime], [],[
#         echo "gomp library is required for this EF5"
#         exit -1])
#fi;
LIBS=-lz -ltiff -lgeotiff
AC_CHECK_HEADERS(xtiffio.h)
AC_CHECK_HEADERS(geotiff/xtiffio.h)
AC_CHECK_HEADERS(libgeotiff/xtiffio.h)
if test x"$ac_cv_header_geotiff_xtiffio_h" = x"yes"; then
        CPPFLAGS="-I/usr/include/geotiff $CPPFLAGS";
else
        if test x"$ac_cv_header_libgeotiff_xtiffio_h" = x"yes"; then
                CPPFLAGS="-I/usr/include/libgeotiff $CPPFLAGS";
        else
                if test x"$ac_cv_header_xtiffio_h" != x"yes"; then
                        AC_MSG_ERROR([xtiffio.h header not found]);
                fi;
        fi;
fi;
AM_INIT_AUTOMAKE([1.9.6 -Wall -Werror no-define foreign])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
