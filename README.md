# Track-Evolution

Track-Evolution is the collection of three numerical approaches to solve the track function evolution in QCD at NLO, which are used in Ref. [1, 2]. In the following, we will describe the usage of moment, Fourier series and Legendre wavelet methods.

### Moment

Moment method solves the track evolution equation in moment space first and then give a polynomial approximation based on moment information. In this case, the moment evolution should be calculated very accurately in order to reduce the error in x-space. 



The moment method evolution is written in the Wolfram language. The LO and NLO evolution in moments are accessible with the functions`evoLO`and `evoNLO ` in the notebook _4QCD_solv_RGEs_RK_release.nb_. These functions use Runge-Kutta fourth-order method to solve the differential equation.



We provide at most 28 moments approximation and the corresponding analytic expression for the kernels are given in the directory _EvolutionEquations4Moments_.

### Fourier Series

### Legendre Wavelet



##### References

[1]  H. Chen, M. Jaarsma, Y. Li, I. Moult, W. Waalewijn, H. X. Zhu, _Collinear Parton Dynamics Beyond DGLAP_, 2210.xxxxx.

[2]  H. Chen, M. Jaarsma, Y. Li, I. Moult, W. Waalewijn, H. X. Zhu, _Multi-Collinear Splitting Kernels for Track Function Evolution_, 2210.xxxxx.
