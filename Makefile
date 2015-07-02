# Master makefile for fastcap and related programs and documentation

SRC = src
FASTH = $(SRC)/fasthenry
ZBUF = $(SRC)/zbuf
MISC = $(SRC)/misc

default:
	@echo Specify what to make:
	@echo " fasthenry - inductance calculation program"
	@echo " zbuf - geometry postscript file program"
	@echo " misc - utility programs"
	@echo " all - all of the above"
	@echo " clean - remove object files"

fasthenry:
	cd $(FASTH) ; $(MAKE) fasthenry

zbuf: 
	cd $(ZBUF) ; $(MAKE) zbuf

misc:
	cd $(MISC) ; $(MAKE) misc

all: fasthenry zbuf misc

clean:
	cd $(FASTH) ; $(MAKE) clean
	cd $(ZBUF) ; $(MAKE) clean
	cd $(MISC) ; $(MAKE) clean

