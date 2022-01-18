# SPM_Project
This is the repository for my SPM project.
For this project I brought 4 parallel implementations of the Recursive Least Square (RLS) online training algorithm for Echo State Networks (ESN). Of these, one uses standard instruments from c++ programming language, while the other three use a library for parallel computation known as FastFlow.  


copied from par0.cpp:
// WHERE TO LOOK TO GET TO THE GIST FAST  
// most interesting flow is:  
// par0.cpp::main-> metrics_computation.hpp::compute_statistics -> metrics_computation.hpp::compute_average_time -> parallel_functs.h::par_train  
// then, the function par_train in file parallel_functs.h distinguishes the various implementations and does the training  

## submodules:
Fisrst of all, the project needs [Fastflow](https://github.com/fastflow/fastflow) library to work.
Then, the project needs [Eigen](https://gitlab.com/libeigen/eigen) and [Spectra](https://github.com/yixuan/spectra/) libraries to run. They are needed for linear algebra that is in turn necessary to generate the contractive recurrent matrix of the reservoir in the ESN

## how to run the project
First of all, one needs to set up the submodules list above.  
Once the dependencies are ready, you can edit the file par0.cpp to edit the parameters you find on the first lines of function main - namely scale of the problem, number of trials per training and maximum parallelism degree - uncomment the parts of the code that you want to run (running all the experiments altogheter takes a lot of time) and then compile by typing on terminal  
make par0  
then run it with  
./par0  
The results are dumped on apposite text files.  

## structure of the project
- The ESN folder contains code for the ESN class 
- The linear_algebra folder contains code that is mainly used for generation of the ESN member objects, It is also needed for statistics-computation and to read dataset
- The [matplotlib-cpp](https://github.com/lava/matplotlib-cpp) submodule was used to do plots from the cpp code. It is not needed anymore but code that used to run it is commented
- The numerical_results_on_report folder conains numerical data that was used to generate the plots you can find on the report 
- The par_utils folder contains the gist of the project: inside you can find the functions for parallel and sequential training togheter with the classes they use 
- BTCUSDT-1m-data.csv contains the dataset. (it reports prices of bitcon in OHLC format)
- par0.cpp is the file with the main function.
- plots_and_stats.ipynb is the python notebook that used to generate the final plots and stats reported on the project. The plots inside were made with data in the numerical_results_on_report folder.
- you can also find the report in report.pdf 

## other stuff
To check that vectorization is correct, have to compile with -fopt-inf-vec-missed and grep the loops in the par_utils/basic_calculations.hpp and tasks.h files. Indeed, these files contain the basic calculations whose vectorization is particularily interesting :

g++ -fdiagnostics-color=always -g par0.cpp linear_algebra/linear_algebra.cpp ESN/ESN.cpp par_utils/utimer.cpp -pthread -lpython3.8 -o par0 -pipe -O3 -I/usr/include/python3.8 -I ./spectra/include -I ./eigen -I ./matplotlib-cpp/ -I./linear_algebra -I./ESN -I./par_utils -I./fastflow -fopt-info-vec-missed 2> >(grep par_utils/basic_calculations.hpp)

g++ -fdiagnostics-color=always -g par0.cpp linear_algebra/linear_algebra.cpp ESN/ESN.cpp par_utils/utimer.cpp -pthread -lpython3.8 -o par0 -pipe -O3 -I/usr/include/python3.8 -I ./spectra/include -I ./eigen -I ./matplotlib-cpp/ -I./linear_algebra -I./ESN -I./par_utils -I./fastflow -fopt-info-vec-missed 2> >(grep par_utils/tasks.h)
