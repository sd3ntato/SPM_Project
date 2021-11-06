main: main.cpp linear_algebra/linear_algebra.cpp linear_algebra/linear_algebra.h ESN/ESN.cpp ESN/ESN.h
	g++ -fdiagnostics-color=always -g \
	main.cpp linear_algebra/linear_algebra.cpp linear_algebra/linear_algebra.h ESN/ESN.cpp ESN/ESN.h \
	-lpython3.8 -o main -pipe -O3 -fopenmp  -ggdb3 \
	-I/usr/include/python3.8 -I ./spectra/include -I ./eigen -I ./matplotlib-cpp/ -I./linear_algebra -I./ESN

seq1: seq1.cpp linear_algebra/linear_algebra.cpp linear_algebra/linear_algebra.h ESN/ESN.cpp ESN/ESN.h
	g++ -fdiagnostics-color=always -g \
	seq1.cpp linear_algebra/linear_algebra.cpp linear_algebra/linear_algebra.h ESN/ESN.cpp ESN/ESN.h \
	-lpython3.8 -o seq1 -pipe -O3 -fopenmp -Wall -ggdb3 \
	-I/usr/include/python3.8 -I ./spectra/include -I ./eigen -I ./matplotlib-cpp/ -I./linear_algebra -I./ESN


clean:
	rm seq1 main