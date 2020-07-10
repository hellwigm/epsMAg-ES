## epsMAg-ES
 Matlab code of the epsMAg-ES for constrained optimization problems

1. ### REFERENCE:  
   
  Hellwig, Michael; Beyer, Hans-Georg, "A Matrix Adaptation Evolution Strategy for Constrained Real-Parameter Optimization", Proceedings of IEEE Conference on Evolutionary Computation (CEC 2018), IEEE Xplore, https://doi.org/10.1109/CEC.2018.8477950, 2018.

2. ### NOTICE:  

  The epsMAg-ES code is made available for reproduction of reported results and testing. The use of the code is permitted subject to the condition of properly acknowledging this source (https://github.com/hellwigm/epsMAg-ES/) as well as citing the relevant papers.

3. ### The epsMAg-ES:  

  The epsMAg-ES represents a novel Evolution Strategy for constrained optimization that combines the the recently suggested Matrix Adaptation Evolution Strategy (MA-ES) with successful constraint handling techniques from the field of Differential Evolution. Being applied to the benchmark problems specified for the CEC 2018 competition on constrained single objective real-parameter optimization, the algorithm is able to find feasible solutions on more than 80% of the benchmark problems with high accuracy. 

4. ### Notice:

  This repository includes the benchmarking functions corresponding to the <b>"CEC 2017/2018 competition on constrained single objective real-parameter optimization"</b> that are openly available in the Githup repository <href>https://github.com/P-N-Suganthan/CEC2017</href>

5. ### Content:

  The pure epsMag-ES algorithm consists of the following modules:

* __epsMAgES.m__ - Main component of the epsMAg-Es algorithm (based on the MA-ES)
* __eps_sort.m__ - Sorting routine for ranking candidate solutions w.r.t. the epsilon-level ordering (setting epsilon to zero results in a lexicographic ordering)
* __eps_rank.m__ - Subroutine of __eps_sort__
* __keep_range.m__ - Box-constraint handling method (reflection into the box)
* __gradientMutation.m__ - Implementation of the gradient-based repair step

  The additional files correspond to the benchmark suite the algorithm was first tested on. The files are included within this repository to provide a demonstration in a running environment. NOTICE: The file __Main_epsMAgES.m__ includes the standard strategy parameter recommendations in the structure array __input__.
* __Main_epsMAgES.m__ - Executable for running the epsMAg-ES on the correspondinng benchmarking problems
* __build_stats.m__ - Postprocessing routine that builds the statistcs to be reported to the CEC 2017 and 2018 competition
* __CEC2017.m__ - Constrained functions specified for the CEC 2017 and 2018 competitions

6.  ### Description:

  The algortihm has the following inputs:
* problem - Matlab structure array that includes a description of the problem
  > problem.constr_fun_name -- Name of the constrained function to be executed, 
  > problem.lb -- vector of the lower parameter vector bounds, 
  > problem.ub -- vector of the upper parameter vector bounds, 
  > problem.gn -- Number of inequality constraints, 
  > problem.hn -- Number of equality constraints.  
"The latter two being needed for the CEC benchmark specific approach for normalizing the constraint violation."

* input   - Matlab structure array specifying all necessary strategy parameters of the algorithm (population sizes, initial mutation strength, learning rates, etc.)
* CEC_fun_o - CEC benchmark specific; some modifications will be necessary to run the algorithm on your own problems

For details refer to the file __Main_epsMAgES.m__  or to the paper.

  During execution, __epsMAg-ES.m__ repeatedly calls up the subroutines
  > __eps_sort.m__, 
  > __eps_rank.m__,  
  > __keep_range.m__,  
  > __gradientMutation.m__, 

that do not need to be individually configured.

  The epsMAg-ES produces the following outputs:
* out - Array of CEC benchmark specific content. Can be omitted in a different context.
* global_best - Structure array containing information of the best candidate solution observed during the algorithm run.
* dyn - Cell array providing information of strategy specific dynamics logged during the run. Capturing these data might be ommitted to reduce execution time.

7. ### License: 
  
  MIT License

Copyright (c) [2018] [Michael Hellwig]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
