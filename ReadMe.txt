MLSGA framework
Copyright (C) 2019  Przemyslaw A. Grudniewski and Adam J. Sobey

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

*******************************
The provided code is based on the scientific research of P.A. Grudniewski and A.J. Sobey [1-3]

*******************************

The following program acts as a complex framework for the Genetic Algorithms (GA) testing and simulation. It allows to utilise a range of genetic algorithms in combination with MLSGA and cMLSGA, or separately. Currently implemented are NSGAII [4], MOEAD [5], MOEADPSF [6], MOEADMSF [6], MOEADM2M [7], MTS [8], BCE [9], HEIA [10] and IBEA [11]. Furthermore, a wide range of benchmarking problems is implemented in order to allow GAs testing. Those includes unconstrained and constrained static problems such as ZDT [12], DTLZ [13], MOP [7], IMB [14], UF [15], WFG [16],  and CF [15], and dynamic cases such as UDF [17], JY [18], FDA [19] and CDF [20]. The different GAs and problems can be operated in loops and all data is stored in accessible format in a single folder. Furthermore, the hyper-parameter tuning can be performed effortlessly due to utilisation of separated loops. The obtained data for each run is presented as the best set of non-dominated solutions in the single file (separately for each algorithm, problem and set of parameters). Furthermore, if chosen the search pattern and the best solutions can be visualised in a form of graphs, videos and heat maps.

The current version of code can compiled only under the MS Visual Studio environment and is compatible with MS Windows OS only. However, the versions Linux version and an unbounded (from VS) code is planned in the future. The detailed instructions of installation and utilisation are provided in Instructions.txt file.

****************Contact***************
You have any questions regarding the code, please contact
Przemyslaw A. Grudniewski	at	pag1c18@soton.ac.uk	or
Adam J. Sobey			at	ajs502@soton.ac.uk


****************Third party content***************
All implemented GAs, such as NSGAII [4], MOEAD [5], MOEADPSF [6], MOEADMSF [6], MOEADM2M [7], MTS [8], BCE [9], HEIA [10] and IBEA [11] are based on the scientific journals articles and corresponding literature can be found in the references. The utilised code for those routines has been based on the provided open-source code or recreated from literature. The respective credits and licensing are given in the corresponding source files.

LibSVM used for the classification is taken as is from https://www.csie.ntu.edu.tw/~cjlin/libsvm/ and used under given copyrights. Those are provided in respective source and header files.

The contour plot generation is based on the Kernel Density Estimation by Tim Nugent (C) 2014 (https://github.com/timnugent/kernel-density)

****************References***************
[1] A. J. Sobey and P. A. Grudniewski, “Re-inspiring the genetic algorithm with multi-level selection theory: Multi-level selection genetic algorithm”, Bioinspiration and Biomimetics, vol. 13, no. 5, pp. 1–13, 2018.

[2] P. A. Grudniewski and A. J. Sobey, “Behaviour of Multi-Level Selection Genetic Algorithm (MLSGA) using different individual-level selection mechanisms”, Swarm Evol. Comput., vol. 44, no. September 2018, pp. 852–862, 2018.

[3] P. A. Grudniewski and A. J. Sobey, "cMLSGA: co-evolutionary Multi-Level Selection Genetic Algorithm", 2019

[4] K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, “A fast and elitist multiobjective genetic algorithm: NSGA-II,” IEEE Trans. Evol. Comput., vol. 6, no. 2, pp. 182–197, 2002.

[5] Q. Zhang and H. Li, “MOEA/D: A Multiobjective Evolutionary Algorithm Based on Decomposition,” IEEE Trans. Evol. Comput., vol. 11, no. 6, pp. 712–731, 2007.

[6] S. Jiang, S. Yang, Y. Wang, and X. Liu, “Scalarizing Functions in Decomposition-Based Multiobjective Evolutionary Algorithms,” IEEE Trans. Evol. Comput., vol. 22, no. 2, pp. 296–313, 2018.

[7] H. Liu, F. Gu, and Q. Zhang, “Decomposition of a Multiobjective Optimization Problem into a Number of Simple Multiobjective Subproblems,” IEEE Trans. Evol. Comput., vol. 18, no. 3, pp. 450–455, 2014.

[8] L.-Y. Tseng and C. Chen, “Multiple Trajectory Search for Multiobjective Optimization,” in 2007 IEEE Congress on Evolutionary Computation (CEC 2007), 2007, pp. 3609–3616.

[9] M. Li, S. Yang, and X. Liu, “Pareto or Non-Pareto: Bi-criterion evolution in multiobjective optimization,” IEEE Trans. Evol. Comput., vol. 20, no. 5, pp. 645–665, 2016.

[10] Q. Lin et al., “A Hybrid Evolutionary Immune Algorithm for Multiobjective Optimization Problems,” IEEE Trans. Evol. Comput., vol. 20, no. 5, pp. 711–729, 2016.

[11] E. Zitzler and K. Simon, “Indicator-Based Selection in Multiobjective Search,” in Parallel Problem Solving from Nature - PPSN VIII, 2004, pp. 832–842.


[12] E. Zitzler, K. Deb, and L. Thiele, “Comparison of Multiobjective Evolutionary Algorithms: Empirical Results,” Evol. Comput., vol. 8, no. 2, pp. 173–195, 2000.

[13] K. Deb, L. Thiele, M. Laumanns, and E. Zitzler, “Scalable Test Problems for Evolutionary Multi-Objective Optimization,” 2001.

[14] H. L. Liu, L. Chen, K. Deb, and E. D. Goodman, “Investigating the effect of imbalance between convergence and diversity in evolutionary multiobjective algorithms,” IEEE Trans. Evol. Comput., vol. 21, no. 3, pp. 408–425, 2017.

[15] Q. Zhang, A. Zhou, S. Zhao, P. N. Suganthan, and W. Liu, “Multiobjective optimization Test Instances for the CEC 2009 Special Session and Competition,” 2009.

[16] S. Huband, L. Barone, L. While, and P. Hingston, “A Scalable Multi-objective Test Problem Toolkit,” in Evolutionary Multi-Criterion Optimization, 2005, no. June, pp. 280–295.

[17] S. Biswas, S. Das, P. N. Suganthan, and C. A. Coello Coello, “Evolutionary Multiobjective Optimization in Dynamic Environments: A Set of Novel Benchmark Functions,” in 2014 IEEE Congress on Evolutionary Computation (CEC 2014), 2014, vol. 1, pp. 3192–3199.

[18] S. Jiang and S. Yang, “Evolutionary Dynamic Multiobjective Optimization: Benchmarks and Algorithm Comparisons,” IEEE Trans. Cybern., vol. 47, no. 1, pp. 198–211, Jan. 2017.

[19] M. Farina, K. Deb, and P. Amato, “Dynamic multiobjective optimization problems: Test cases, approximations, and applications,” IEEE Trans. Evol. Comput., vol. 8, no. 5, pp. 425–442, 2004.

[20] P. A. Grudniewski and A. J. Sobey, “Genetic Algorithm performances on constrained dynamic problem types”, 2019.
