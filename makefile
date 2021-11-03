main: main.cpp linear_algebra.cpp linear_algebra.h ESN.cpp ESN.h
	rm main;\
	g++ -fdiagnostics-color=always -g *.cpp -I/usr/include/python3.8 -lpython3.8 -o main -pipe -O3 -fopenmp 

clean:
	rm main