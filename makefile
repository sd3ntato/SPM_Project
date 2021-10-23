main: main.cpp linear_algebra.cpp linear_algebra.h ESN.cpp ESN.h
	rm main;\
	g++ -fdiagnostics-color=always -g *.cpp -o main

clean:
	rm main