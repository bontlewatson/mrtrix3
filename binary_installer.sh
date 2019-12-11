#!/bin/bash -ex
export PATH=/bin:/sbin:/usr/bin:/usr/sbin
PLD=$(pwd)
BRANCH=dev
THREADS=$(sysctl -n hw.ncpu)
TMPDIR=$(mktemp -d)
echo ${TMPDIR}
cd ${TMPDIR}
PREFIX=${TMPDIR}/mrtrix3

if [ -f "${PLD}/mrtrix3-build-deps.tar.gz" ]; then
  tar xfz ${PLD}/mrtrix3-build-deps.tar.gz
else
  # Grab latest versions from git
  EIGEN_VERSION=$(git ls-remote --tags https://github.com/eigenteam/eigen-git-mirror.git | awk '{print $2}' | grep -v '\^{}$' | grep -v alpha | grep -v beta | grep -v '-' | sort -V | tail -1)
  EIGEN_VERSION=${EIGEN_VERSION#*/*/}
  echo "Using eigen ${EIGEN_VERSION}"

  TIFF_VERSION=$(git ls-remote --tags https://gitlab.com/libtiff/libtiff.git | awk '{print $2}' | grep -v '\^{}$' | grep -v alpha | grep -v beta | grep v.* | sort -V | tail -1)
  TIFF_VERSION=${TIFF_VERSION#*/*/v}
  echo "Using tiff ${TIFF_VERSION}"
  
  PNG_VERSION=$(git ls-remote --tags https://github.com/glennrp/libpng.git | awk '{print $2}' | grep -v '\^{}$' | grep -v alpha | grep -v beta | grep v.* | sort -V | tail -1)
  PNG_VERSION=${PNG_VERSION#*/*/v}
  echo "Using png ${PNG_VERSION}"
  
  FFTW_VERSION=$(git ls-remote --tags https://github.com/FFTW/fftw3.git | awk '{print $2}' | grep -v '\^{}$' | grep -v alpha | grep -v beta | grep 'fftw-3' | sort -V | tail -1)
  FFTW_VERSION=${FFTW_VERSION#*/*/fftw-}
  echo "Using fftw ${FFTW_VERSION}"
  
  QT_VERSION=$(git ls-remote --tags https://github.com/qt/qt5.git | awk '{print $2}' | grep -v '\^{}$' | grep -v alpha | grep -v beta | grep -v '-' | grep '5.9' | sort -V | tail -1)
  QT_VERSION=${QT_VERSION#*/*/v}
  echo "Using qt ${QT_VERSION}"

  # EIGEN
  SECONDS=0
  curl -s -O -L http://bitbucket.org/eigen/eigen/get/${EIGEN_VERSION}.tar.bz2
  tar xf ${EIGEN_VERSION}.tar.bz2
  mkdir -p ${PREFIX}/include/eigen3
  cp -R eigen*/Eigen ${PREFIX}/include
  cp -R eigen*/unsupported ${PREFIX}/include
  EIGEN_SECONDS=${SECONDS}

  # TIFF
  SECONDS=0
  curl -s -O -L http://download.osgeo.org/libtiff/tiff-${TIFF_VERSION}.tar.gz
  tar xf tiff-${TIFF_VERSION}.tar.gz
  cd tiff-${TIFF_VERSION}
  ./configure -q -prefix ${PREFIX} --enable-shared=NO --without-x
  make install > /dev/null
  cd ..
  TIFF_SECONDS=${SECONDS}

  # PNG
  SECONDS=0
  curl -s -O -L https://ftp.osuosl.org/pub/blfs/conglomeration/libpng/libpng-${PNG_VERSION}.tar.xz
  tar xf libpng-${PNG_VERSION}.tar.xz
  cd libpng-${PNG_VERSION}
  ./configure -q -prefix ${PREFIX} --enable-shared=NO
  make install > /dev/null
  cd ..
  PNG_SECONDS=${SECONDS}

  # FFTW
  SECONDS=0
  curl -s -O http://www.fftw.org/fftw-${FFTW_VERSION}.tar.gz
  tar xf fftw-${FFTW_VERSION}.tar.gz
  cd fftw-${FFTW_VERSION}
  ./configure -q -prefix ${PREFIX} --disable-doc --disable-fortran --disable-debug --enable-threads --disable-dependency-tracking --enable-sse2 --enable-avx
  make install > /dev/null
  cd ..
  FFTW_SECONDS=${SECONDS}

  # QT5 BASE
  SECONDS=0
  curl -s -O http://ftp1.nluug.nl/languages/qt/archive/qt/${QT_VERSION%.*}/${QT_VERSION}/submodules/qtbase-opensource-src-${QT_VERSION}.tar.xz
  tar xf qtbase-opensource-src-${QT_VERSION}.tar.xz
  cd qtbase-opensource-src-${QT_VERSION}
  ./configure -opensource -confirm-license -release -no-dbus -no-openssl -no-harfbuzz -no-freetype  -no-cups -no-framework -nomake examples -prefix ${PREFIX}
  make -j ${THREADS}
  make install
  cd ..
  QTBASE_SECONDS=${SECONDS}

  # QT5 SVG
  SECONDS=0
  curl -s -O http://ftp1.nluug.nl/languages/qt/archive/qt/${QT_VERSION%.*}/${QT_VERSION}/submodules/qtsvg-opensource-src-${QT_VERSION}.tar.xz
  tar xf qtsvg-opensource-src-${QT_VERSION}.tar.xz
  cd qtsvg-opensource-src-${QT_VERSION}
  ${PREFIX}/bin/qmake
  make -j ${THREADS}
  make install
  cd ..
  QTSVG_SECONDS=${SECONDS}

  tar cpfz ${PLD}/mrtrix3-build-deps.tar.gz mrtrix3
fi

# MRTRIX3
SECONDS=0
git clone https://github.com/MRtrix3/mrtrix3.git mrtrix3-src -b ${BRANCH}
cd mrtrix3-src
CFLAGS="-I${PREFIX}/include" LINKFLAGS="-L${PREFIX}/lib" TIFF_LINKFLAGS="-llzma -ltiff" PATH=${PREFIX}/bin:${PATH} ./configure
NUMBER_OF_PROCESSORS=${THREADS} ./build
MRTRIX_VERSION=$(cat lib/mrtrix3/_version.py | awk '{print $3}' | tr -d '"')
cd ..
MRTRIX_SECONDS=${SECONDS}


mv ${PREFIX} ${PREFIX}_dep
mkdir -p ${PREFIX}/mrtrix3
cp -R ${PREFIX}-src/bin    ${PREFIX}/mrtrix3
cp -R ${PREFIX}-src/lib    ${PREFIX}/mrtrix3
cp -R ${PREFIX}-src/share  ${PREFIX}/mrtrix3
cp -R ${PREFIX}-src/matlab ${PREFIX}/mrtrix3
cp ${PREFIX}-src/set_path  ${PREFIX}/mrtrix3
cp ${PREFIX}_dep/lib/libQt5{Core,Gui,OpenGL,PrintSupport,Svg,Widgets}${QT_POSTFIX}.*.dylib ${PREFIX}/mrtrix3/lib
mkdir -p ${PREFIX}/mrtrix3/bin/plugins/{platforms,imageformats}
cp ${PREFIX}_dep/plugins/platforms/libqcocoa${QT_POSTFIX}.dylib  ${PREFIX}/mrtrix3/bin/plugins/platforms
cp ${PREFIX}_dep/plugins/imageformats/libqsvg${QT_POSTFIX}.dylib ${PREFIX}/mrtrix3/bin/plugins/imageformats

cp -R ${PLD}/MRView.app ${PREFIX} 
mkdir -p ${PREFIX}/MRView.app/Contents/MacOS/
mv ${PREFIX}/mrtrix3/bin/mrview ${PREFIX}/MRView.app/Contents/MacOS/
cp ${PLD}/mrview ${PREFIX}/mrtrix3/bin

cp -R ${PLD}/SHView.app ${PREFIX}     
mkdir -p ${PREFIX}/SHView.app/Contents/MacOS/
mv ${PREFIX}/mrtrix3/bin/shview ${PREFIX}/SHView.app/Contents/MacOS/
cp ${PLD}/shview ${PREFIX}/mrtrix3/bin

tar cfz ${PLD}/mrtrix3.tar.gz ${PREFIX}

TOTAL_SECONDS=$((EIGEN_SECONDS + TIFF_SECONDS + PNG_SECONDS + FFTW_SECONDS + QTBASE_SECONDS + QTSVG_SECONDS + MRTRIX_SECONDS))
echo "eigen ${EIGEN_VERSION}: ${EIGEN_SECONDS} s"
echo "tiff ${TIFF_VERSION}: ${TIFF_SECONDS} s"
echo "png ${PNG_VERSION}: ${PNG_SECONDS} s"
echo "fftw ${FFTW_VERSION}: ${FFTW_SECONDS} s"
echo "qtbase ${QT_VERSION}: ${QTBASE_SECONDS} s"
echo "qtsvg ${QT_VERSION}: ${QTSVG_SECONDS} s"
echo "mrtrix ${MRTRIX_VERSION}: ${MRTRIX_SECONDS} s"
echo "total : ${TOTAL_SECONDS} s"
rm -rf ${TMPDIR}
