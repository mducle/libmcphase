# Makefile for all McPhase

bindir = .
include $(bindir)/src/Makefile.common

cfdir = $(bindir)/cf1ion_module
icdir = $(bindir)/ic1ion_module
phdir = $(bindir)/phonon_module
bfkdir = $(bindir)/bfk_src
mcpdir = $(bindir)/src
vecdir = $(bindir)/src/vector
funcdir = $(bindir)/src/functions
calcdir = $(bindir)/src/calc

all: vector functions mcphase phonon

vector: 
	cd $(vecdir) && $(MAKE) 

functions:
	cd $(funcdir) && $(MAKE)

calc:
	cd $(calcdir) && $(MAKE)

cfield: vector
	cd $(cfdir) && $(MAKE)

ic1ion: vector
	cd $(icdir) && $(MAKE)

bfk   : 
	cd $(bfkdir) && $(MAKE)

phonon   : 
	cd $(phdir) && $(MAKE)

mcphase: vector ic1ion cfield
	cd $(mcpdir) && $(MAKE)

clean: 
	cd $(vecdir) && $(MAKE) clean
	cd $(funcdir) && $(MAKE) clean
	cd $(cfdir) && $(MAKE) clean
	cd $(icdir) && $(MAKE) cleanall
	cd $(phdir) && $(MAKE) cleanall
	cd $(mcpdir) && $(MAKE) clean
	cd $(phdir) && $(MAKE) clean
	rm -vf $(bindir)/addj.exe $(bindir)/charges.exe $(bindir)/coq2jjj.exe \
		$(bindir)/mcdispit.exe $(bindir)/singleion.exe $(bindir)/cfield.exe \
		$(bindir)/cond.exe $(bindir)/jjj2j.exe $(bindir)/mcphasit.exe $(bindir)/spins.exe \
                $(bindir)/chrgplt.exe $(bindir)/pointc.exe $(bindir)/spinsfromq.exe \
                $(bindir)/mcdiff.exe $(bindir)/cf1ion_module/cfield.dll \
                $(bindir)/ic1ion.exe $(bindir)/icf1ion.exe $(bindir)/so1ion.exe \
                $(bindir)/ic1ion_module/ic1ion.dll \
                $(bindir)/fediff.exe $(bindir)/mf2fe.exe \
		$(bindir)/formfactor.exe $(bindir)/radwavfunc.exe
	rm -vf $(bindir)/addj $(bindir)/charges $(bindir)/coq2jjj \
		$(bindir)/mcdispit $(bindir)/singleion $(bindir)/cfield \
		$(bindir)/cond $(bindir)/jjj2j $(bindir)/mcphasit $(bindir)/spins \
                $(bindir)/chrgplt $(bindir)/pointc $(bindir)/spinsfromq \
                $(bindir)/mcdiff $(bindir)/cf1ion_module/cfield.so \
                $(bindir)/ic1ion $(bindir)/icf1ion $(bindir)/so1ion \
                $(bindir)/ic1ion_module/ic1ion.so $(bindir)/ic1ion_module/icf1ion.so \
                $(bindir)/fediff $(bindir)/mf2fe \
                $(bindir)/formfactor $(bindir)/radwavfunc
