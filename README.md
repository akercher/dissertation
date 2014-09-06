dissertation
============

Source code written for dissertation.

athena4.2
  - exact_mhd.c : exact_solver for ideal MHD.
    Install GNU scientific library (http://www.gnu.org/software/gsl/).  If necessary, edit lines
      GSLINC = -I/usr/include
      GSLLIB = -L/usr/lib/ -lgsl -lgslcblas
    of Makeoptions.in if local paths are different. 

latex
  - dissertation and defense source files.