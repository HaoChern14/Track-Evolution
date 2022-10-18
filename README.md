# Track-Evolution

Track-Evolution is the collection of three numerical approaches to solve the track function evolution in QCD at NLO, which are used in Ref. [1, 2]. In the following, we will describe the usage of moment, Fourier series and Legendre wavelet methods.

### Moment

Moment method solves the track evolution equation in moment space first and then give a polynomial approximation based on moment information. In this case, the moment evolution should be calculated very accurately in order to reduce the error in x-space. 



The moment method evolution is written in the Wolfram language. The LO and NLO evolution in moments are accessible with the functions`evoLO` and `evoNLO` in the notebook `4QCD_solv_RGEs_RK_release.nb`. These functions use Runge-Kutta fourth-order method to solve the differential equation.



We provide at most 28 moments approximation and the corresponding analytic expression for the kernels are given in the directory `EvolutionEquations4Moments`.

### Fourier Series

Fourier series method solves the evolution equation in the Fourier coefficient space and becomes a finite dimensional ODE after truncation to finite Fourier modes. The numerical evolution kernels for Fourier method are stored in the [datasets hosted on Zenodo](https://zenodo.org/record/7219729#.Y063VOxBz-Q), where we provide at most 100-mode kernels at LO and 40-mode kernels at NLO. The Julia script `NLO_evolution_fourier.jl` makes use of the Julia package _DifferentialEquations_ to solve the ODE.



For easier use, we create [a docker image](https://hub.docker.com/r/haochern/qcd-track-evolution-fourier) which packs the Fourier method evolution script and the evolution kernels.  The [Fouier docker image](https://hub.docker.com/r/haochern/qcd-track-evolution-fourier) should be used under the folder `Fourier/NLO_evolution_fourier` in which the initial conditions are provided and final results are stored. In this part, we will assume `./ ` to be `Fourier/NLO_evolution_fourier/`.

##### Provide initial conditions and energy scales

- The file `./data/initial_scale.txt`  is the energy scale (in the unit GeV) of the initial conditions. 

- In the `./data/target_scales.txt`, we can input any reasonable number of wanted energy scales to see the evolution results. 

- `./data/Tq_initial_coeff_re.txt` and `./data/Tq_initial_coeff_im.txt` are the real and imaginary part of quark track function Fourier coefficients lists 
  $$
  \int_{0}^{1} e^{-2\pi i n x} T_q(x)dx, n=1,2,\dots,100\,.
  $$
  The same applies to the gluon case.

##### Solving the differential equations

0. Download [the docker image](https://hub.docker.com/r/haochern/qcd-track-evolution-fourier):

   ```
   docker pull haochern/qcd-track-evolution-fourier:1.0
   ```
   
1. Execute
   
   ```
   cd data;
   docker run -v ${PWD}:/home/data haochern/qcd-track-evolution-fourier:1.0
   ```
   or ```./run.sh``` in the direcory `Fourier/NLO_evolution_fourier/`.

##### Plot the results
The Fourier docker image offers two kinds of output -- both the Fourier coefficients results and discretized curves. The Mathematica notebook `plot.nb` is an example to manipulate the output data.

In the folder `data/output/`, `NLO_fouier_coeff_mu_xxx.txt` is the Fourier coefficients output at energy scale `xxx`, where the first 100 complex numbers are the quark track function Fourier coefficients and the next 100 complex numbers are for the gluon case. `Tq_value_mu_xxx.txt` and `Tg_value_mu_xxx.txt`are the values of track functions evaluated at the discretized points given in `x_axis.txt`.

   

### Legendre Wavelet



##### References

[1]  H. Chen, M. Jaarsma, Y. Li, I. Moult, W. Waalewijn, H. X. Zhu, _Collinear Parton Dynamics Beyond DGLAP_, 2210.xxxxx.

[2]  H. Chen, M. Jaarsma, Y. Li, I. Moult, W. Waalewijn, H. X. Zhu, _Multi-Collinear Splitting Kernels for Track Function Evolution_, 2210.xxxxx.
