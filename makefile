main: main.cpp linear_algebra.cpp linear_algebra.h ESN.cpp ESN.h
	rm main;\
	g++ -fdiagnostics-color=always -g *.cpp -lpython3.8 -o main -pipe -O3 -fopenmp -Wall -ggdb3 \
	-I/usr/include/python3.8 -I ./spectra/include -I ./eigen -I ./matplotlib-cpp/


clean:
	rm main