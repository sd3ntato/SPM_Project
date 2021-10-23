main: main.cpp funct.cpp funct.h ESN.cpp ESN.h
	rm main;\
	g++ -fdiagnostics-color=always -g *.cpp -o main

clean:
	rm main