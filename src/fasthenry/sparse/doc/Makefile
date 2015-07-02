# Makefile for Sparse documentation
#

PRINTER = -Pprado
# PRINTER = -Plps40

default		: usersguide
all		: usersguide install advert

usersguide	: spDoc.ms
	tbl spDoc.ms | eqn | ptroff -ms $(PRINTER)

install		: install.ms
	ptroff -ms $(PRINTER) install.ms

advert		: advert.ms
	tbl advert.ms | ptroff -ms $(PRINTER)

all-lpr		: usersguide.lpr install.lpr

usersguide.lpr	: spDoc
spDoc		: spDoc.ms
	tbl spDoc.ms | neqn | nroff -ms | colcrt - > spDoc

install.lpr	: install.ms
	nroff -ms install.ms | colcrt - > install.lpr

