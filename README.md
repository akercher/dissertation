dissertation
============

* This dissertation includes extensions to version 4.2 of Athena, a code for astrophysical MHD.  
  The official version of Athena can be downloaded from
  https://trac.princeton.edu/Athena/wiki/AthenaDocsDownLd. 

athena4.2
  - exact_mhd.c : exact_solver for ideal MHD.
    Install GNU scientific library (http://www.gnu.org/software/gsl/).  If necessary, edit lines
      GSLINC = -I/usr/include
      GSLLIB = -L/usr/lib/ -lgsl -lgslcblas
    of Makeoptions.in if local paths are different. 

latex
  - dissertation and defense source files.