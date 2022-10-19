# Track-Evolution

Track-Evolution is the collection of three numerical approaches to solve the track function evolution in QCD at NLO, which are used in Ref. [1, 2]. In the following, we will describe the usage of moment, Fourier series and Legendre wavelet methods. The existing example of initial conditions are extracted from Pythia at 100 GeV [3]. 

## Moment

Moment method solves the track evolution equation in moment space first and then give a polynomial approximation based on moment information. In this case, the moment evolution should be calculated very accurately in order to reduce the error in x-space. 



The moment method evolution is written in the Wolfram language. The LO and NLO evolution in moments are accessible with the functions`evoLO` and `evoNLO` in the notebook `4QCD_solv_RGEs_RK.nb`, which use Runge-Kutta fourth-order method to solve the differential equation. More usage details are explained in the notebook.



We provide at most 28 moments approximation and the corresponding analytic expression for the kernels are given in the directory `EvolutionEquations4Moments`.

## Fourier Series

Fourier series method solves the evolution equation in the Fourier coefficient space and becomes a finite dimensional ODE after truncation to finite Fourier modes. The numerical evolution kernels for Fourier method are stored in the [datasets hosted on Zenodo](https://zenodo.org/record/7219729#.Y063VOxBz-Q), where we provide at most 100-mode kernels at LO and 40-mode kernels at NLO. The Julia script `NLO_evolution_fourier.jl` makes use of the Julia package _DifferentialEquations_ to solve the ODE.



For easier use, we create [a docker image](https://hub.docker.com/r/haochern/qcd-track-evolution-fourier) which packs the Fourier method evolution script and the evolution kernels.  The [Fouier docker image](https://hub.docker.com/r/haochern/qcd-track-evolution-fourier) should be used under the folder `Fourier/NLO_evolution_fourier` in which the initial conditions are provided and final results are stored. In this part, we will assume `./ ` to be `Fourier/NLO_evolution_fourier/`.

##### Provide initial conditions and energy scales

- The file `./data/initial_scale.txt`  is the energy scale (in the unit GeV) of the initial conditions. 

- In the `./data/target_scales.txt`, we can input any reasonable number of wanted energy scales to see the evolution results. 

- `./data/Tq_initial_coeff_re.txt` and `./data/Tq_initial_coeff_im.txt` are the real and imaginary part of quark track function Fourier coefficients lists 
  $$b^q_n=\int_{0}^{1} e^{-2\pi i n x} T_q(x)dx,\; n=1,2,\dots,100\,.$$
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
   or 
   
   ```./run.sh``` in the direcory `Fourier/NLO_evolution_fourier/`.

##### Plot the results
The Fourier docker image offers two kinds of outputs -- both the Fourier coefficients results and discretized curves. The Mathematica notebook `plot.nb` is an example to manipulate the output data.

In the folder `data/output/`, `Tq_value_mu_xxx.txt` and `Tg_value_mu_xxx.txt`are the values of track functions evaluated at the discretized points given in `x_axis.txt` and energy scale `xxx`.
 `NLO_fouier_coeff_mu_xxx.txt` is the Fourier coefficients output at energy scale `xxx`, where the first 100 complex numbers $\{b^q_1, \dots, b^q_{100}\}$ are the quark track function Fourier coefficients and the next 100 complex numbers $\{b^g_1, \dots, b^g_{100}\}$ are the gluon case. To revover the x-space track functions $T_{i=q,g}(x)$, we can simply use
 $$T_i(x) = 1 + 2\,\mathrm{Re}\sum_{i=1}^{100} b^i_n e^{2\pi i n x}\,.$$
   

## Legendre Wavelet
Our wavelet method approximates the track functions with piecewise polynomials. Specifically, we divide the range $(0,1)$ into 16 intervals, on each of which we use the quadratic polynomial approximation. This corresponds to set parameters $K=5$ and $M=2$ in Sec. 7.1.2 in the Ref. [2]. The numerical kernels can be found on [the same Zenodo website](https://zenodo.org/record/7219729#.Y063VOxBz-Q) as that for the Fourier method.

We also provide [the wavelet docker image](https://hub.docker.com/r/haochern/qcd-track-evolution-wavelet) which should be used with the folder `Wavelet/NLO_evolution_wavelet`. Therefore, we assume `./` = `Wavelet/NLO_evolution_wavelet/` in this part. Its usage is similar to the Fourier docker image.

##### Provide initial conditions and energy scales

- The file `./data/initial_scale.txt`  is the energy scale (in the unit GeV) of the initial conditions. 

- In the `./data/target_scales.txt`, we can input any reasonable number of wanted energy scales to see the evolution results. 

- `./data/Tq_initial_coeff.txt`are the wavelet coefficients for quark track function.
  The same applies to the gluon case. An example of generating the wavelet coefficient lists from given functions is shown in the Mathematica notebook `intitial_coefficients_generation_example.nb`.

##### Solving the differential equations

0. Download [the docker image](https://hub.docker.com/r/haochern/qcd-track-evolution-wavelet):

   ```
   docker pull haochern/qcd-track-evolution-wavelet:1.0
   ```
   
1. Execute
   
   ```
   cd data;
   docker run -v ${PWD}:/home/data haochern/qcd-track-evolution-wavelet:1.0
   ```
   or 
   
   ```./run.sh``` in the direcory `Wavelet/NLO_evolution_wavelet/`.

##### Plot the results
The wavelet docker image outputs are `data/output/NLO_wavelet_coeff_mu_xxx.txt`, which stores the quark and gluon wavelet coefficients at energy scale `xxx`. An example of plotting the track functions with wavelet coefficients is given in the Mathematica notebook `plot.nb`.

#### References

[1]  H. Chen, M. Jaarsma, Y. Li, I. Moult, W. Waalewijn, H. X. Zhu, _Collinear Parton Dynamics Beyond DGLAP_, 2210.xxxxx.

[2]  H. Chen, M. Jaarsma, Y. Li, I. Moult, W. Waalewijn, H. X. Zhu, _Multi-Collinear Splitting Kernels for Track Function Evolution_, 2210.xxxxx.

[3]  H.-M. Chang, M. Procura, J. Thaler, and W. J. Waalewijn, _Calculating track thrust with track functions_, [Phys. Rev. D 88, 034030](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.88.034030), [[1303.6637]](https://arxiv.org/abs/1303.6637).
