all : TestLinearInterpolator

TestLinearInterpolator : TestLinearInterpolator.cpp LinearInterpolator.h
	clang++ -o $@ $<
