SWIG = swig
SWIG_OPTS = -c++ -python
CC = g++
OPTS = -fPIC -funroll-loops -O3 -msse3 -fopenmp

all: swig1 swig2 swig3 link

swig1: M2Pwrap.i
	$(SWIG) $(SWIG_OPTS) $?

swig2: P2Pwrap.i
	$(SWIG) $(SWIG_OPTS) $?

swig3: calculateMultipoles.i
	$(SWIG) $(SWIG_OPTS) $?

link:  
	$(CC) $(OPTS) -c M2Pwrap.cpp M2Pwrap_wrap.cxx -I/usr/include/python2.7 -I/usr/lib/python2.7
	$(CC) $(OPTS) -fPIC -shared $? -o _M2Pwrap.so M2Pwrap.o M2Pwrap_wrap.o
	$(CC) $(OPTS) -c P2Pwrap.cpp P2Pwrap_wrap.cxx -I/usr/include/python2.7 -I/usr/lib/python2.7
	$(CC) $(OPTS) -fPIC -shared $? -o _P2Pwrap.so P2Pwrap.o P2Pwrap_wrap.o
	$(CC) $(OPTS) -c calculateMultipoles.cpp calculateMultipoles_wrap.cxx -I/usr/include/python2.7 -I/usr/lib/python2.7
	$(CC) $(OPTS) -fPIC -shared $? -o _calculateMultipoles.so calculateMultipoles.o calculateMultipoles_wrap.o

clean:
	rm -f _M2Pwrap.so M2Pwrap.o M2Pwrap_wrap.o M2Pwrap_wrap.cxx M2Pwrap.py *.pyc
	rm -f _P2Pwrap.so P2Pwrap.o P2Pwrap_wrap.o P2Pwrap_wrap.cxx P2Pwrap.py *.pyc
	rm -f _calculateMultipoles.so calculateMultipoles.o calculateMultipoles_wrap.o calculateMultipoles_wrap.cxx calculateMultipoles.py *.pyc

.cxx.o:
	$(CC) $(OPTS) -c $< -o $@

.cpp.o:
	$(CC) $(OPTS) -c $< -o $@
