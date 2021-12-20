prova: prova.cpp \
  linear_algebra/linear_algebra.cpp linear_algebra/linear_algebra.h \
	ESN/ESN.cpp ESN/ESN.h \
	par_utils/utimer.cpp \
	par_utils/tasks.cpp par_utils/tasks.h \
	par_utils/pool.cpp par_utils/pool.h \
	par_utils/prot_queue.h
	g++ -fdiagnostics-color=always -g \
	prova.cpp \
	linear_algebra/linear_algebra.cpp linear_algebra/linear_algebra.h \
	ESN/ESN.cpp ESN/ESN.h \
	par_utils/utimer.cpp \
	par_utils/tasks.cpp par_utils/tasks.h \
	par_utils/pool.cpp par_utils/pool.h \
	par_utils/prot_queue.h \
	-pthread \
	-lpython3.8 -o prova -pipe -O3 -fopenmp  -ggdb3 \
	-I/usr/include/python3.8 -I ./spectra/include -I ./eigen -I ./matplotlib-cpp/ -I./linear_algebra -I./ESN -I./par_utils

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
	rm seq1 main prova

par0: par0.cpp \
  linear_algebra/linear_algebra.cpp linear_algebra/linear_algebra.h \
	ESN/ESN.cpp ESN/ESN.h \
	par_utils/utimer.cpp \
	par_utils/ff_Pool.hpp \
	par_utils/tasks.h \
	par_utils/pool.h \
	par_utils/parallel_functs.h \
	par_utils/seq_functs.h \
	par_utils/utils.h \
	par_utils/prot_queue.h
	g++ -fdiagnostics-color=always -g \
	par0.cpp \
	linear_algebra/linear_algebra.cpp \
	ESN/ESN.cpp \
	par_utils/utimer.cpp \
	-pthread \
	-lpython3.8 -o par0 -pipe -O3 \
	-I/usr/include/python3.8 -I ./spectra/include -I ./eigen -I ./matplotlib-cpp/ -I./linear_algebra -I./ESN -I./par_utils -I./fastflow

prova1: prova1.cpp 
	g++ prova1.cpp -O3 -o prova1 -I ./par_utils/ -I ./linear_algebra/ -I ./eigen/ -I ./spectra/include/ -I ./fastflow/ -pthread

