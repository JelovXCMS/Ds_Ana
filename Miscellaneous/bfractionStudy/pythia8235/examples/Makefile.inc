# PYTHIA configuration file.
# Generated on Wed Jun 20 10:10:40 EDT 2018 with the user supplied options:
# --with-root=/Users/liangzhang/root

# Install directory prefixes.
PREFIX_BIN=/Users/liangzhang/Desktop/pythia8235/bin
PREFIX_INCLUDE=/Users/liangzhang/Desktop/pythia8235/include
PREFIX_LIB=/Users/liangzhang/Desktop/pythia8235/lib
PREFIX_SHARE=/Users/liangzhang/Desktop/pythia8235/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
ENABLE_SHARED=false
CXX=g++
CXX_COMMON=-O2  -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-dynamiclib
CXX_SONAME=-Wl,-dylib_install_name,@rpath/
LIB_SUFFIX=.dylib

# EVTGEN configuration.
EVTGEN_USE=false
EVTGEN_BIN=
EVTGEN_INCLUDE=./
EVTGEN_LIB=./

# FASTJET3 configuration.
FASTJET3_USE=false
FASTJET3_BIN=
FASTJET3_INCLUDE=./
FASTJET3_LIB=./

# HEPMC2 configuration.
HEPMC2_USE=false
HEPMC2_BIN=
HEPMC2_INCLUDE=./
HEPMC2_LIB=./

# HEPMC3 configuration.
HEPMC3_USE=false
HEPMC3_BIN=
HEPMC3_INCLUDE=./
HEPMC3_LIB=./

# LHAPDF5 configuration.
LHAPDF5_USE=false
LHAPDF5_BIN=
LHAPDF5_INCLUDE=./
LHAPDF5_LIB=./

# LHAPDF6 configuration.
LHAPDF6_USE=false
LHAPDF6_BIN=
LHAPDF6_INCLUDE=./
LHAPDF6_LIB=./

# POWHEG configuration.
POWHEG_USE=false
POWHEG_BIN=
POWHEG_INCLUDE=./
POWHEG_LIB=./

# PROMC configuration.
PROMC_USE=false
PROMC_BIN=
PROMC_INCLUDE=./
PROMC_LIB=./

# ROOT configuration.
ROOT_USE=true
ROOT_BIN=/Users/liangzhang/root/bin/
ROOT_INCLUDE=/Users/liangzhang/root/include
ROOT_LIB=/Users/liangzhang/root/lib

# GZIP configuration.
GZIP_USE=false
GZIP_BIN=
GZIP_INCLUDE=./
GZIP_LIB=./

# BOOST configuration.
BOOST_USE=false
BOOST_BIN=
BOOST_INCLUDE=./
BOOST_LIB=./

# PYTHON configuration.
PYTHON_USE=false
PYTHON_BIN=
PYTHON_INCLUDE=./
PYTHON_LIB=./
