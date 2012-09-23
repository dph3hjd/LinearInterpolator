SCI_INCLUDES:=-I/Applications/Scientific/install/include
SCI_LIBS:=-L/Applications/Scientific/install/lib -lTTextToData -lTFieldSelection
ROOT_CONFIG:=/Applications/Scientific/install/bin/root-config #$(shell which root-config)
ROOT_CFLAGS=$(shell $(ROOT_CONFIG) --cflags)
ROOT_LIBS=$(shell $(ROOT_CONFIG) --libs --glibs --auxlibs)


all : TestLinearInterpolator

TestLinearInterpolator : TestLinearInterpolator.cpp LinearInterpolator.h
	clang++ -o $@ $(SCI_INCLUDES) $(SCI_LIBS) $(ROOT_CFLAGS) $(ROOT_LIBS) $<
