# Master makefile for fastcap and related programs and documentation

SRC = src
FASTH = $(SRC)/fasthenry
ZBUF = $(SRC)/zbuf

default:
	@echo Specify what to make:
	@echo " fasthenry - inductance calculation program"
	@echo " zbuf - geometry postscript file program"
	@echo " all - all of the above"
	@echo " clean - remove object files"

fasthenry:
	cd $(FASTH) ; $(MAKE) fasthenry

zbuf: 
	cd $(ZBUF) ; $(MAKE) zbuf

all: fasthenry zbuf

clean:
	cd $(FASTH) ; $(MAKE) clean
	cd $(ZBUF) ; $(MAKE) clean
