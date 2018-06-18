### Package Makefile ###

PKGDIR=~/Documents/synthsys/yeast-stress-network/modelling/MIdecoding
PKGNAME=MIdecoding

all: install

.PHONY: package
package: check install doc

.PHONY: install
install: build
	./Rpkgscript install ${PKGDIR} ${PKGNAME}

.PHONY: check
check: build
	./Rpkgscript check ${PKGDIR} ${PKGNAME}

.PHONY: build
build:
	./Rpkgscript build ${PKGDIR} ${PKGNAME}

.PHONY: doc
doc:
	./Rpkgscript doc ${PKGDIR} ${PKGNAME}
