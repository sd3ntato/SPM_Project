# SPM_Project
This is the repository for my SPM project.
For this project I brought 4 parallel implementations of the Recursive Least Square (RLS) online training algorithm for Echo State Networks (ESN). Of these, one uses standard instruments from c++ programming language, while the other three use a library for parallel computation known as FastFlow.
## submodules:
Fisrst of all, the project needs [Fastflow](https://github.com/fastflow/fastflow) library to work.
Then, the project needs [Eigen](https://gitlab.com/libeigen/eigen) and [Spectra](https://github.com/yixuan/spectra/) libraries to run. They are needed for linear algebra that is in turn necessary to generate the contractive recurrent matrix of the reservoir in the ESN
## how to run the project
Once the dependencies are ready, you can edit the file par0.cpp to edit the parameters you find on the first lines of function main - namely scale of the problem, number of trials per training and maximum parallelism degree - uncomment the parts of the code that you want to run (running all the experiments altogheter takes a lot of time) and then type on terminal

make par0

to compile the executable and then run it with 

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