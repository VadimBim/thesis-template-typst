#import "thesis_template.typ": *
#import "common/titlepage.typ": *
#import "thesis_typ/disclaimer.typ": *
#import "thesis_typ/acknowledgement.typ": *
#import "thesis_typ/abstract_en.typ": *
#import "common/metadata.typ": *

//Packages
#import "@preview/lovelace:0.2.0": *
#import "@preview/physica:0.9.2": *

#titlepage(
  title: titleEnglish,
  degree: degree,
  program: program,
  supervisor: supervisor,
  advisors: advisors,
  author: author,
)

#disclaimer(
  title: titleEnglish,
  degree: degree,
  author: author,
)

#acknowledgement()

#abstract_en()

#show: project.with(
  title: titleEnglish,
  degree: degree,
  program: program,
  supervisor: supervisor,
  advisors: advisors,
  author: author,
)

= Bayesian Optimization <chapter-bo>

#emph("Bayesian Optimization") (BO) is a class of machine-learning-based algorithms that are designed to optimize expensive, black-box functions $f: A arrow.r RR$ (from this point forward, we will refer to these functions as #emph("objective functions")). Usually, these algorithms are efficient when the objective has the follwing properties @Frazier2018-dq:

- The input set is a hypter-rectangle $A = {vb(x) in RR^d | a_i lt.eq x_i lt.eq b_i}$ where $d$ is not large (typically, $d lt.eq 20$).
- The objective function is continuous.
- $f$ is "expensive to evaluate". This is especially true in the cases where objective is a result of a computationally demanding simulation, where each run can take hours.
- $f$ is a "black box" function, which means that we lack any knowledge about the structure of the function (linearity, concavity, periodicity, etc.). If we knew some of the properties, we could use more efficient algorithms that leverege this information.
- When we evaluate our objective, we observe only $f(vb(x))$ and no-first or second-derivatives.
- $f$ can be noisy, which is very relevant in the case where objective is obtained from an physical experiment.
- Also, BO is suited for cases where we want to find a global, rather than a local extremum.

The algorithm is composed of two main components: a statistical #emph("surrogate model") of the objective and an #emph("acqusition function"). The common choice for surrogate is to use Gaussian Processes (shortly GP) regression (@f_sim_gp) to model the function that we want to optimize. However, it can be anything that can make a prediction with the associated uncertainty about it.

$ f(vb(x)) tilde cal(G) cal(P) (mu(vb(x)), k(vb(x), vb(x)')) $ <f_sim_gp>

Informally, we can read @f_sim_gp as follow: the objective function follows a normal distribution over possible functions, determined by #emph("mean function") $mu(vb(x))$ and #emph("covariance function") $k(vb(x), vb(x)')$ (also known as the #emph("kernel")). 

The acqusition function $alpha(vb(x))$ is a heuristic that "tells" us where is the most promising point to evaluate our objective at the next iteration. The most promising point is given by the extremum of $alpha (vb(x))$:

$ vb(x)_n = "argmax" alpha_(n-1) (vb(x)) $

These functions are "cheaper" to optimize relative to the objective, so that traditional optimization techniques are employed. There are a lot of acqusition functions based on different assumptions. For example, a widely used acqusition function is Expected Improvement $upright(E)upright(I)(vb(x))$, as the name suggests this function returns the expected improvement relative to the best observed value so far.. 

Basic Bayesian Optimization is a sequential algorithm, here the pseudo-code is shown (adapted from @Frazier2018-dq):

#show: setup-lovelace

#algorithm(
  caption: [Basic pseudo-code for Bayesian Optimization],
  pseudocode(
    [Place a Gaussian process prior on $f$],
    [Observe $f$ at $n$ initial points. $"Data" := cal(D) = {(vb(x)_1, y_1), ..., (vb(x)_n, y_n)}$],
    [*while* $n lt.eq N$ *do*], ind,
      [Update the posterior probability distribution on $f$ using $cal(D)$],
      [$vb(x)_(n + 1) = "argmax" alpha_(n)(vb(x))$],
      [Observe $y_(n+1) = f(vb(x)_(n+1))$],
      [Add $(vb(x)_(n+1), y_(n+1))$ to $cal(D)$], ded,
    [*end while*],
    [*return* best $y$]
  )
)

First $n$ points are usually generated using Sobol sequences @Sobol1967-nw to evenly sample from input space. Here $N$ is the number of available evaluations. 

In the next sections we disucss the components of BO in more details.

== Gaussian Processes <chapter-gp>

The formal defintion is as follow @Rasmussen2005-ou:

#rect(
  width: 100%,
  radius: 10%,
  stroke: 0.5pt,
  fill: rgb("#E6F9FF"),
)[#emph("A Gaussian process is a collection of random variables, any finite number of which have a joint Gaussian distribution")]

More specifically, if $f$ is sampled from a GP, then a finite number function values $vb(f) := [f(vb(x)_1), f(vb(x)_2), ..., f(vb(x)_n)]^T$ follows a multivariate normal distribution:

$ vb(f) &tilde cal(N) (vb(mu), vb(K)) = (2 pi)^(-n slash 2) |vb(K)|^(-1 slash 2) exp ( - 1/2 (vb(f) - vb(mu))^T vb(K)^(-1) (vb(f) - vb(mu))) $

, where $vb(mu)$ is the mean vector and $vb(K)$ is the $n times n$ covariance matrix given by the kernel $K_(i j) = k(vb(x)_i , vb(x)_j)$. Kernels have the property that points closer in input space are more strongly correlated: $ norm(vb(x) - vb(x)') < norm(vb(x) - vb(x)'') arrow.r.double k(vb(x), vb(x)') > k(vb(x), vb(x)'') $

If we want to make a new observation $f(vb(x)_(n+1)) := f_(n+1)$, by defintion of GP, it will come from the same probability distribution as $vb(f)$. Thus $P(f_(n+1) | vb(f))$ is obtained from marginalization (@marginalization) #footnote[All figures in @chapter-gp were generated using  http://www.infinitecuriosity.org/vizgp/.] of the underlying joint distribution $P([f_1, ... , f_n, f_(n+1)]^T)$, which is also a multivariate normal distribution. In conlusion. it can be shown @Rasmussen2005-ou that: // TODO put the derivation in the appendix.

$ &P(f_(n+1) | vb(f)) tilde cal(N) (mu_(n+1), sigma_(n+1)^(space 2)), "where" \
  &mu_(n+1) = k(x_(n+1), vb(x)) vb(K)^(-1) (vb(f) - vb(mu)) + mu(x_(n+1)) \
  &sigma_(n+1)^(space 2) = k(x_(n+1), x_(n+1)) - k(x_(n+1), vb(x)) vb(K)^(-1) k(vb(x), x_(n+1)) \
  &k(x_(n+1), vb(x)) := [k(x_(n+1), x_1), ..., k(x_(n+1), x_n)] = k(vb(x), x_(n+1))^T $ <posterior-dist>

#figure(
  image("figures/marginalization.png", width: 85%),
  caption: [The marginal distribution (#text(orange)[orange]) $P(f(2) | 3 "observations")$.],
) <marginalization>

In @marginalization we have an unkown 1D function, where we made 3 evaluations. Based on available data, we infered the probability distribution of $f(2)$.

The covariance matrix should be symmetric and positive, which limits the number of possible functions that $k$ can take. If we want to model smooth functions, the usual choice for the kernel is the Squared Exponential kernel (SE):

$ k(vb(x), vb(x')) &= a exp(- norm(vb(x) - vb(x'))^2), \ 
  norm(vb(x) - vb(x'))^2 &:= sum_(i=1)^d (x_i - x'_i) / (l_i^2), $ <SE_kernel>

$vb(theta) := [a, l_1, ..., l_d]$ are known as the #emph("hyperparameters") of the kernel, and they are responsible of the form of the sampled function. It is interesting to note that these hyperparameters are interpretable (see @samples_se). For example, if we take into consideration case $d=1$, we have only 2 hyperparameters $a$ and $l_1$. In this case, $a$ will influence the amplitude of the sampled functions, and $l_1$ the lengthscale (how fast our functions vary). 

#figure(
table(
  stroke: none,
  columns: 2,
  gutter: 1,
  align: horizon,
  inset: (x:5pt, y:20pt),
  image("figures/samples_1_1.png", width: 105%),
  image("figures/samples_1,4_0,2.png", width: 100%),
  ),
  kind: image,
  caption: [4 sampled functions from GP with a SE kernel (bottom) and associated slices through the kernel $k(x_i, dot.c)$ (top) as a function of the second argument. On the left side $[a, l_1^2] = [1, 2]$ and on the right $[a, l_1^2] = [1.4, 0.08]$. ┈ are means, #box(rect(fill: rgb("#A2C6C6"), height: 8pt, width: 8pt, radius: 1pt), baseline: 5%) represents $plus.minus sigma$ and #box(rect(fill: rgb("#CADEDE"), height: 8pt, width: 8pt, radius: 1pt), baseline: 5%) $plus.minus 2 sigma$ confidence bands]
) <samples_se>

As we mentioned before, kernels are responsible for the structure of the modeled functions. Morevoer, we can combine different kernels via addition and multiplication to obtain functions with specific characteristics (see Chapter 2 of @Duvenaud2014).

#figure(
  image("figures/se_times_linear.png", width: 80%),
  caption: [Locally periodic 1D samples drawn from a GP with a  SE $times$ Per kernel. Where Per is a periodic kernel $:= a exp (-2 (sin^2(pi |x - x'| slash p ))/l^2)$]
)

For the mean function, the most common choice is a #emph("const.") value. It is possible to construct a mean that captures a specific trend of the function:

$ mu(vb(x)) = italic("const.") + sum_(i=1)^p beta_i psi_i (vb(x)), $

where $psi_i (vb(x))$ are known parametric functions, usually low-order polynomials. However, in the case of BO, we don't have paticular knowledge about the objective.

The modeled functions are very sensitive to the chosen hyperparameters, and the natural question is: How do we choose these hyperparameters such that sampled functions describe the best observed data? We will mention only Bayesian model selection, as the rest is out of scope of the current study (interested reader can consult Chapter 5 of @Rasmussen2005-ou for details). By Bayes rule, observed function values and hyperparameters are related as follow #footnote[Here we assumed that model is fixed by specifying form of the kernel and mean function, otherwise we should have written: $P(vb(theta), cal(M) | vb(f))$.]:

$ P(vb(theta) | vb(f)) = (P(vb(f) | vb(theta)) P(vb(theta))) / P(vb(f)) $

$P(vb(theta) | vb(f))$ is known as the #emph("posterior"), $P(vb(f) | vb(theta))$ is the #emph("likelihood") and $P(vb(f))$ is a normalization #emph("const.") named #emph("marginal likelihood") (a.k.a. evidence):

$ P(vb(f)) = integral P(vb(f) | vb(theta)) P(vb(theta)) dd(vb(theta)) $

 Most of the times, this integral is not analytically tractable and analytical approximation or Monte Carlo methods are used. As mentioned before, we want to chhose hyperparameters that maximize the probability to see observed values :

$ hat(vb(theta)) = arg max P(vb(theta) | vb(f)) = arg max P(vb(f) | vb(theta)) P(vb(theta)), $ <MAP>

where we used the fact that $P(vb(f))$ doesn't depend on $vb(theta)$. In @MAP we are estimating hyperparameters by #emph("Maximum a posteriori") (MAP) estimate. An approximation can be obtained by assuming that $P(vb(theta))$ has a constant density over $vb(theta)$, such #emph("Maximum likelihood estimation") (MLE) (@MLE) is obtained. One can see that MAP is more computationally intensive, however it is a good choice when MLE gives unreasonable hyperparameters.

$ hat(vb(theta)) = arg max P(vb(f) | vb(theta)) $ <MLE>

Since GP assumes that function values are generated from a multivariate normal distribution, and taking the natural log we have:

$ ln P(vb(f) | vb(theta)) = -n/2 ln 2 pi - 1/2 ln |vb(K)| - 1/2 (vb(f) - vb(mu))^T vb(K)^(-1) (vb(f) - vb(mu)) $ <logP>

An important detail is that if our observations are noisy $vb(y) := vb(f) + vb(epsilon)$, we can model it by adding a diagonal matrix to the covariance matrix. In the case in which we assume that all obervations have the same variance of the noise (#emph("homoscedastic") noise), we simply add $sigma^2 I$ to $vb(K)$. When noise is different for each observation (#emph("heterorscedastic") noise), associated variances are added to the diagonal. Moreover, we can add correlations between noises to off-diagonal elements if known.

// TODO: derivation eq: 5.9 din @Rasmussen2005

Because $vb(theta)$ in @logP is nested inside the correlation matrix, deriving an analytical formula is not feasible. Therefore, we need to employ numerical optimization techinques to maximize $ln P(vb(f) | vb(theta))$.

== Acqusition functions <chapter-af>

In this section we will discuss in detail the most commonly used acqusition function, namely Expected Improvement. We will limit our analysis to 1D case and also will mention other functions used for #emph("standard") BO problems (check chapter 5 of @Frazier2018-dq for an introduction in Exotic BO).

Suppose that we performed a number of evaluations and $f(x^*)$ is the best value observed so far. Now we have one additional evaluation to perform, and we can perform it anywhere. We define improvement as $I (x) := max(f(x) - f(x^*), 0)$. While we would like to choose $x$ such that this improvement is large, $f(x)$ is unkown until after the evaluation. What we can do, however, is to take the expected value of this improvement and choose $x$ to maximize it. We define the expected improvement:

$ "EI" (x) := E [I (x)], $ <EI_definition>

where $"E"[dot.c]$ indicates the expectation taken under the posterior distribution . Expected improvement can be written in closed form (see @A_EI):

$ "EI" (x) &= (mu(x) - f(x^*)) (1 - Phi(z_0)) + sigma phi(z_0), \
  z_0 &= (f(x^*) - mu(x))/sigma(x) $ 

Another widely used acqusition function is #emph("Probability of Improvement"):

$ "PI" (x) := P(f(x) > f(x^*)). $ <PI_definition>

which can be simply evaluated with the same change of variables as in @A_EI, giving us :

$ "PI" (x) = Phi ((mu(x) - f(x^*)) / sigma(x)) $

The difference between PI and EI is that PI only considers the probability of improvement and not the expected magnitude of the new evaluation. The simplest acqusition function is #emph("Upper Confidence Bound"):

$ "UCB" (x) := mu(x) + beta sigma(x), $

where $beta$ controls the exploration/exploitaion tradeoff. When $beta$ is small, the algorithm will explore regions with high $mu(x)$, on the contrary, when $beta$ is large, BO rewards the exploration of currently uncharted areas (high $sigma(x)$). A much simpler approach is #emph("Thompson sampling"), where we basically draw a sample from the GP and evaluate our function at the maximum of the drawn sample.

#figure(
  image("figures/acq_funcs.jpg", width: 85%),
  caption: [Examples of acqusition functions for standard BO.  \
  #text(0.6em)[ Image source: https://www.borealisai.com/wp-content/uploads/2020/06/T9_2-2.png] ]
)

As we have seen, evaluating an acquisition function requires evaluating an integral over the posterior. However, in higher-dimensional settings, this is usually analytically intractable. An alternative is to use Monte-Carlo sampling to approximate the integrals @Monte-Catlo-botorch. For example, we can approximate EI as follow :

$ "EI" (vb(X)) approx 1/N sum_(i=1)^N max_(j = 1, ..., q) I(xi_(i j)), $ <EI_monte-carlo>

where $vb(X) = [vb(x)_1, ... , vb(x)_q]$ is discretized input space and $vb(xi)_i$ is a realization of $P([f(vb(x)_1), ..., f(vb(x)_q)])$. Usually $q > N$.

#figure(
  image("figures/expected_improvement_mc.png", width: 85%),
  caption: [3 Monte Carlo realizations $vb(xi)_i$, represented by dotted lines, given #text(rgb("#99C1F1"))[observed data]. #text(red)[Red] horizontal lines represents argument of the maximum value of improvement for a monte carlo sample]
)

= Laser Wakefield Acceleration (LWFA)

LWFA, for the first time proposed by Tajima and Dawson in 1979 @Tajima1979, is a promising technique for obtaining high-energy electron bunches on acceleration distances much smaller compared to the conventional accelerators (see @RF_vs_plasmcavity). The main idea is to use high electric fields obtained when an ultra-short, high-power laser pulse goes through an #emph("underdense") (more details in the following sections) plasma. Through the effect of the #emph("ponderomotive force"), electrons are expelled, while ions that are much heavier than electrons ($m_i >> m_e$) can be considered immobile. This charge imbalance creates a longitudinal plasma wave that follows the driver (laser). In the bubble-like structure behind the driver, electrons can be injected using different methods (self-injection, optical injection, ionization injection, etc.), each having its own advantages and disadvantages.

#figure(
  image("figures/RF_vs_plasmcavity.png", width: 85%),
  caption: [Left: Readiofrequency cavity. Right: Laser plasma Wakefield. The paser pulse in #text(yellow)[yellow] propagates from left to right, the iso-electric density os hown in #text(blue)[blue] and the lectron bunch in red #text(red)[red]. Image credits to @Malka2016]
) <RF_vs_plasmcavity>

In the next sections, we will define plasma and how we can describe it. Then, we will introduce plasma waves and the ponderomotive force. We will end this chapter with the description of the #emph("bubble") regime of LWFA, which is suited for acceleration of electrons up to relativistic energies.

== Plasma 

#rect(
  width: 100%,
  radius: 10%,
  stroke: 0.5pt,
  fill: rgb("#E6F9FF"),
)[#emph("Plasma is a quasi-neutral ionized gas that can support collective phenomena.")]

#emph("Quasi-neutrality") means that charge densities of positive and negative particles are approximately equal $n_(+) tilde.eq n_(-)$. Usually, most negatively charged particles are electrons, and positively charged ones are ions. The density of particles and temperature $T$ determine the type of plasma (see @plasma-types). In the LWFA scenario, we are dealing with "cold" and "dense" plasma. Cold means that the range for the kinetic energy of electrons is $E_c tilde.eq [10^2, 10^3] "eV"$. Dense means that we are in the density range $n_e = 10^(18) - 10^(24) "cm"^(-3)$. In @plasma-types, $N_D$ is the number of electrons in the #emph("Debye sphere") (see @plasma-lt-scale).

#pagebreak()

#figure(
  image("figures/Plasma-types.png", width: 85%),
  caption: [Example of plasma types in the density-temperature plane. Image credits to @Gibon2020]
) <plasma-types>

The most accurate and feasible description of plasma dynamics is captured by kinetic theory, which describes the evolution of #emph("distribution function") $f_s (vb(r), vb(p), t)$ in the phase space. The distribution function shows the density of particles of species $s$ in the voxel centered at point $(vb(r), vb(p))$ in phase space at time $t$. Physical observables can be obtained by integrating moments of $f_s$. For example, if we integrate this object over momentum coordinate, we can obtain the average density for each species $n_s (vb(r), t)$ and mean velocity $vb(u)_s (vb(r), t)$:

#grid(columns: (1fr, 1fr), math.equation(block: true, numbering: none, $ n_s (vb(r), t) = integral f_s dd(vb(p)^3) $), [$ vb(u)_s (vb(r), t) = 1/n_s integral vb(v) f_s dd(vb(p)^3) $ <currents>])

For a collisionless #footnote([If we want to include collisions in our model, we are inserting a #emph("collision integral") in the RHS of @Vlasov-eq]) plasma, the continuity equation in phase space is:

$ pdv(f_s, t) + vb(v) dprod grad_(vb(r)) f_s + q_s (vb(E) + vb(v) cprod vb(B)) dprod grad_(vb(p)) f_s = 0, $ <Vlasov-eq>

where $vb(v) = vb(p)/(m_s gamma_s) = (c vb(p))/sqrt(m_s^2 c^2 + vb(p)^2)$ is the velocity of particles in relativistic regime. The above equation is known as the Vlasov equation. Together with the Maxwell equations, it forms a #emph("self-consistent") system of equations that describes plasma dynamics. Where sources are given by:

#grid(columns: (1fr, 1fr), math.equation(block: true, numbering: none, $ rho = sum_s q_s n_s $), [$ vb(j) = sum_s q_s n_s vb(u)_s $ <sources>])

Where $q_s$ is the charge of species $s$. However, the Vlasov-Maxwell system is impossible to solve analytically except for a limited number of simple systems. In our case, we will employ a simpler model derived from kinetic theory (see section 2.2.2 @Macchi2013), the so-called #emph("fluid model") of plasma. this model assumes that each species with density $n_s$ and velocity $vb(u)_s$ behaves in a fluid-like manner. We will consider ion fluid stationary, and for the electrons we have:

$ &pdv(n_e, t) + div (n_e vb(u)_e) = 0 \
  &n_e m_e dv(vb(u)_e, t) = - n_e e (vb(E) + vb(u)_e cprod vb(B)) - grad P_e \
  &"Maxwell equations" \
  &dv(, t) (P_e n_e^(-gamma_e)) = 0 $ <fluid_model>

The last equation in @fluid_model is known as the closure condition, making our equations complete. $gamma_e = (2 + N)/N$ is the specific heat ratio, where $N$ is the number of degrees of freedom for the electron.

== Waves in plasma

Plasma, like any other fluid, can transfer energy through waves. We will derive the #emph("dispersion relation") of plasma waves by following the approach in section 3.1 from @Gibon2020. 

Let's consider @fluid_model in 1D case. We will note with $Z$ the number of protons an ion has. Also, we will suppose that fluid velocities are small, making the term $ prop vb(u)_e cprod vb(B) tilde.eq vb(0)$. Because our system is constrained to 1D, electrons have only 1 degree of freedom $N=1 => gamma_e = 3$. With all of this, one can write:

$ pdv(n_e, t) + pdv((n_e u_e), x) &= 0 \ 
  n_e (pdv(u_e, t) + u_e pdv(u_e, x)) &= -e/m n_e E - 1/m pdv(P_e, x) \
  dv(, t) ( P_e / (n_e^3) ) &= 0  \ 
  pdv(E, x) = e/epsilon_0 (Z n_i - n_e) &underbracket(=, Z n_i tilde.eq n_0) e/epsilon_0 (n_0 - n_e) $ <fluid_model_1D>

Now, we have 4 non-linear differential equations. To linearize them, we will consider that our system exhibits small perturbations and will ignore higher-order terms in these perturbations:

$ n_e := n_0 + n_1 \
  u_e := u_1 \
  P_e := P_0 + P_1 \
  E := E_1. $ <1D_perturbations>

Taking this into account, @fluid_model_1D becomes:

$ pdv(n_1, t) + n_0 pdv(u_1, x) &= 0 \
  n_0 pdv(u_1, t) &= - e/m_e n_0 E_1 - 1/m_e pdv(P_1, x) \ 
  pdv(E_1, x) &= - e / epsilon_0 n_1 \
  P_1 &= 3 k_B T_e n_1 $ <fluid_model_1D_linearized>

For pressure, we assumed isothermal background electrons $P_0 = k_B T_e n_0$ (check @linearization_pressure for details). If we take time derivative of the first equation in @fluid_model_1D_linearized, write $pdv(u_1, t, s: \/)$ using second equation and eliminate $pdv(E_1, x, s: \/)$ and $P_1$ with third and fourth equation, we will get :

$ (pdv(, t, 2) - 3 v_("th")^2 pdv(, x, 2) + omega_p^2) n_1 = 0, $

where $omega_p$ is plasma frequency as defined in @plasma-lt-scale and $v_("th")^2 := k_B T_e slash m_e$ is the squared of electronic thermal velocity. If we look for plane wave solutions of the form $n_1 = n_0 exp(i(omega t - k x))$, differential operators become: $partial_t arrow.r i omega$ and $partial_x arrow.r -i k$, giving us in the end the Bohm-Gross dispersion relation:

$ omega^2 = omega_p^2 + 3 v^2_("th") k^2 $ <Bohm-Gross-dispersion>

To derive the dispersion relation for EM transverse waves of the #emph("laser"), we will apply $curl$ to $curl vb(E) = -diff_t vb(B)$ and will use $curl vb(B) = mu_0 epsilon_0 diff_t vb(E) + mu_0 vb(j)$. Where the charge current $vb(j) tilde.eq -n_0 e vb(u)$. Using the vectorial identity $curl curl vb(a) = grad(div vb(a)) - laplacian vb(a)$ and taking into consideration that $div vb(E) = 0$ for pure EM waves with $m_e diff_t vb(u) = -e vb(E) $ (group and phase velocities of the laser are $>> v_("th")$ and we can assume $P_e tilde.eq 0$), we have:

$ (laplacian - 1/c^2 pdv( ,t, 2) - omega_p^2/c^2 ) vb(E) = vb(0) $

With plane wave ansatz, we arrive at the dispersion relation for EM waves:

$ omega_L^2 = c^2 k_L^2 + omega_p^2. $ <tranverse-dispersion>

The dispersion relations give an overview of which propagation modes are permitted.
// TODO Cetz plot with dispersion relation as in @Gibon2020.

Now we can introduce #emph("critical density") by setting $omega_L = omega_p$:

$ n_(e,c) = (m_e epsilon_0)/e^2 omega_L^2 = 1.1/(lambda_L^2[mu m]) 10^(21) "cm"^(-3) $

If $omega_L < omega_p$, the plasma response is fast enough to reflect the incident laser. We call this type of plasma #emph("overdense") because the laser cannot pass through. On the contrary, if $omega_L > omega_p$, the laser can pass, and the plasma is called #emph("underdense"). We now can define phase velocity and plasma velocity for the laser:

#grid(columns: (1fr, 1fr), math.equation(block: true, numbering: none, $ v_phi := omega_L/k_L = c/eta $), $ v_g := pdv(omega_L, k_L) = eta c $)

where $eta = sqrt(1 - omega_p^2/omega_L^2)$ is the refractive index. We can also associate a Lorentz gamma factor as $gamma_L := 1 slash sqrt(1 - (v_g / c)^2) tilde.eq omega_L/omega_p$. One important observation is that the phase velocity of the plasma wave is approximately equal to the group velocity of the laser ($v_phi^p tilde.eq v_g^L$). Indeed,  if we start from @Bohm-Gross-dispersion:

$ v_phi^p := omega/k_p = c sqrt(1 + underbrace((3 v_("th")^2)/c^2, << 1)) tilde.eq c, $

we arrive at the conclusion that it is the same as the group velocity of the laser for an underdense plasma ($omega_p^2 << omega_L^2$).

In the following, we will follow section 3.3 of @Karsch2020 to arrive at an equation that quantitatively describes the excitation of plasma waves. 

We will start by rewriting the second equation from @fluid_model in the relativistic regime ($m_e vb(u)_e arrow.r vb(p) = m_e gamma vb(u)_e$) and will keep $P_e = 0$ because of the same assumption used to derive @tranverse-dispersion. Also, we will write EM fields in terms of potentials $(Phi, vb(A))$ and drop the subindex $#hide("")_e$ for convenience:

$ (diff_t + vb(u) dprod grad) vb(p) = e (grad Phi + diff_t vb(A) - vb(u) cprod curl vb(A) ). $

We can put $e vb(u) cprod curl vb(A)$ to the left and add and subtract $vb(u) cprod (curl vb(p))$. Now, we can use the relation $grad p^2 = 2 [(vb(p) dprod grad) vb(p) + vb(p) (curl vb(p))]$ with $vb(p) = gamma m_e vb(u)$ and $gamma = sqrt(1 + (p/(m_e c))^2)$, to arrive at: 

$ m_e c^2 gamma grad gamma = (vb(u) dprod grad) vb(p) + vb(u) cprod (curl vb(p)) $.

We will consider that the electrons are driven only by the vector potential, which will allow us to set $vb(u) cprod curl(e vb(A) - vb(p)) = 0$. By introducing normalized quantities:

#grid(columns: (1fr, 1fr, 1fr), math.equation(block: true, numbering: none, $ vb(a) := (e vb(A))/(m_e c) $), math.equation(block:true, numbering: none, $ phi := (e Phi)/(m_e c^2) $), math.equation(block: true, numbering: none,$ vb(pi) := vb(p)/(m_e c) $))

,we arrive at:

$ pdv(vb(pi), t) = pdv(vb(a), t) + c grad (phi - gamma). $

The equation states that the dynamics of the electron fluid is dictated by the Coulomb force associated with the charge distribution ($c grad phi$) and the relativistic ponderomotive force ($- c grad gamma$).This equation is also the starting point for exploring the solution of linear and non-linearly driven plasma waves. The regime is given by the parameter $a_0$ #footnote([$a_0$ can be interpreted as the ration of the mean quiver energy to the electron rest mass.]) of the laser ($a_0 > 1 arrow.r$ non-linear, linear otherwise). 

Although analytical treatment of the 3D non-linear regime is not possible, this study will explore it using PIC simulations. In this scenario, the so-called bubble regime can be achieved.

== Bubble regime

Bubble regime is obtained when a spherical cavity free of electrons is formed behind the laser (see @bubble_regime). The laser must be powerful enough ($a_0 gt.tilde 2$) to expel #emph("all") the electrons and to form the ion cavity. To enter this regime, one must fulfill the following matching conditions for the laser intensity, spot size $w_b$, pulse length ($c tau_0 lt.tilde w_b$), and plasma wavelength:

$ w_b tilde.eq r_b = (2 sqrt(a_0))/k_p $

#figure(
  image("figures/bubble_regime.png", width: 60%),
  caption: [Bubble regime observed in a 3D OSIRIS PIC simulation of LWFA acceleration. The laser is propagating towards the right. It can be seen that the bubble is fully developed and that the electron density drops to zero. Created fields exceed the wavebreaking limit, and in the back of the bubble, the transverse self-infection is visible. The energy of the injected electrons is color-coded in black to red. Image adapted from @Karsch2020]
) <bubble_regime>

The radius of this cavity was obtained from phenomenological studies of the 3D PIC simulations and is given by $r_b approx 2 sqrt(a_0) slash k_p$. Also, the fields were estimated for an ideal sphere in this case:

#grid(columns: (1fr, 1fr, 1fr), math.equation(block: true, numbering: none, $ E_z/E_0 = k_p/2 xi $), math.equation(block:true, numbering: none, $ E_r/E_0 = k_p/4 r $), [$ B_theta/E_0 = - k_p /(4 c) r, $ <EM_bubble_field>])

where $E_0 = m_e omega_p c slash e$ is the #emph([cold non-relativistic wavebreaking field]) (plasma fluid velocities are non-relativistic) and $xi := z - v_g t$ is the comoving longitudinal position. This is the value of the electric field at which background electrons can outrun the wake and, therefore, be injected. This process is called longitudinal wavebreaking. For fields greater than $E_0$, the fluid model of plasma @fluid_model is no longer valid. The maximum value of the longitudinal electric field is given by:

$ E_(z,max) = E_0 sqrt(a_0) <=> E_(z,max)["GV"/"m"] tilde.eq 96 sqrt(n_(e, 0) [10^18 "cm"^(-3)]) sqrt(a_0) $

It is important to keep in mind that EM fields from above were obtained with the following main physical assumptions:

1. Matching conditions are perfectly met.
2. The cavity is #strong([completely]) free of electrons.
3. There is no interaction between the electrons bunch and the laser pulse.

As a consequence, the EM fields are overestimated. One can substitute @EM_bubble_field into the equation of motion for the electron, which interacts with EM field $dot(vb(p)) = -e (vb(E) + vb(v) cprod vb(B))$, to derive forces:

#math.equation(block: true, numbering: none, $ &dot(vb(p)) tilde.eq  -(e E_0 k_p)/2 [(r/2 vb(e)_r + xi vb(e)_z) - r/(2c) v_z underbrace(vb(e)_z cprod vb(e)_theta, -vb(e)_r)] =  \
  &= -(m_e omega_p^2) / 2 (1/2(1 + v_z/c) r vb(e)_r + xi vb(e)_z) tilde.eq -(m_e omega_p^2) / 2 (r vb(e)_r + xi vb(e)_z) <=> $)

#grid(columns: (1fr, 1fr), math.equation(block: true, numbering: none, $ F_perp = -(m_e omega_p^2) / 2 r  $), math.equation(block:true, $ F_parallel = -(m_e omega_p^2) / 2 xi $))

We assumed that $v_z >> v_perp$ and $v_z slash c tilde.eq 1$. Transverse focusing of the electrons happens over the whole bubble, while defocusing occurs only on-axis density peak. For $xi in (0, r_b slash 2)$ electrons are decelerating. The distance, relative to the lab frame, over which electrons must propagate before they reach the middle of the bubble is named #emph("dephasing length") $:= L_d >> r_b$. This, with the pump depletion, are the limiting factors in the acceleration process. 

= Betatron Radiation

The accelerated relativistic electrons obtained with LWFA can wiggle strongly. This wiggling causes electrons to radiate energy. The radiation emitted by an electron in the direction of observation $vb(n)$ (assuming that observation of the radiation is made #emph("far") from the electron) is given by (see eq. 14.65 from @Jackson1998-cw):

$ dd(cal(W), 2) / (dd(omega) dd(Omega)) = e^2/(4 pi^2 c) abs(integral_(-oo)^(+oo) dd(t) e^(i omega (t - (vb(n) dprod vb(r)(t))/c )) (vb(n) cprod [(vb(n) - vb(beta)) cprod dot(vb(beta))])/(1 - vb(beta) dprod vb(n))^2)^2. $ <dW_dwdOmgega>

$dd(cal(W), 2) / (dd(omega) dd(Omega))$ is the energy radiated within a spectral band $dd(omega)$ centered on the frequency $omega$ in solid angle $dd(Omega)$ centered on the direction of observation $vb(n)$. $vb(beta)$ is the velocity of the electron normalized to the speed of light $c$, and $vb(r)(t)$ is the position of the electron function of time. Some important conclusions can be derived from the above equation:

1. If $dot(vb(beta)) = 0$, no radiation is observed, which means that acceleration is responsible for the radiation emission.
2. Energy is maximum when $vb(beta) dprod vb(n) arrow.r 1 <=> beta tilde.eq 1$ and $vb(beta) || vb(n)$. Relativistic electrons will radiate orders of magnitude higher than non-relativistic ones. This is a direct consequence of Lorentz transformation. 
3. Because $dot(vb(beta)_(||)) prop vb(F)_(||)/gamma^3$ and $dot(vb(beta)_(perp)) prop vb(F)_(perp)/gamma$, applying transverse force is more efficient for the relativistic electron ($gamma >> 1$) if we want it to radiate stronger.
4. Locally, $e^(i omega (t - vb(n) dprod vb(r)/c)) tilde.eq e^(i omega (1 - beta))$. The integral will give a non-zero result when the integrand (without exponential) oscillates approximatively with the same frequency as the exponential $omega (1 -beta)$. If we define $omega_(e^-)$ the frequency at which $vb(beta)$ varies, for a non-zero result $omega_(e^-) tilde.eq omega (1 - beta)$. The electron will radiate higher at the frequency $omega = omega_(e^-)/(1 - beta) tilde.eq 2 gamma^2 omega_(e^-)$. This is exactly the Doppler upshift. We can obtain X-rays by wiggling a relativistic electron at a frequency far below the X-ray range: $omega_(e^-) tilde.eq omega_X/(2 gamma^2)$.

Suppose the electron follows a simple sinusoidal trajectory $x(z) = x_0 sin(k_u z)$ with a constant velocity $beta$ and the spatial period $lambda_u$.

#figure(
  image("figures/undulator_vs_wiggler.png", width: 55%),
  caption: [Two radiation regimes. The #text(rgb("#C0C3E2"))[shaded] cones represent the direction of the instantaneously emitted radiation with an opening $Delta theta ~ 1 slash gamma$. $Psi$ is the maximum angle between $vb(beta)$ and $vb(e)_z$. Image adapted from @Corde2013 section II. A.]
) <undulator_vs_wiggler>

We can define the parameter $K := gamma Psi$, which tells us in what regime are we (#emph("undulator") for $K << 1$ or #emph("wiggler") $K >> 1$). In the undulator regime, the radiation spectrum is centered at fundamental frequency $omega$, which depends on the angle of observation $theta$. The radiation opening (divergence) in this case is $theta_r = 1/gamma$ 

For the wiggler regime, the spectrum will contain harmonics up to #emph("critical frequency") $omega_c$. The divergence of the emitted radiation in the direction of the oscillation will be $theta_X = K/gamma$ and will remain the same in the direction perpendicular to the oscillation $theta_Y = 1/gamma$.

Now let us consider an ideal electron bunch consisting of $N_e$ electrons. By ideal, we mean that electrons have the same energy and initial momentum. Also, we assume that electrons follow the same reference trajectory $vb(r)(t)$, but with a spatio-temporal translation $(t_j, vb(R)_j)$. In this case, emitted radiation can be written as:

$ dd(cal(W), 2) / (dd(omega) dd(Omega)) = underbrace(abs(sum_(j=1)^N_e e^(i omega (t_j - (vb(n) dprod vb(R)_j)/c)))^2, c(omega)) dprod eval(dd(cal(W), 2) / (dd(omega) dd(Omega)))_("ref"), $

where $eval(dd(cal(W), 2) / (dd(omega) dd(Omega)))_("ref")$ is the radiated energy of the reference trajectory given by @dW_dwdOmgega. $c(omega)$ is the #emph("coherance factor") which depends on the initial distribution of translations $eval({t_j, vb(R)_j})_(j=1)^N_e$. If the sample $eval({t_j, vb(R)_j})_(j=1)^N_e$ follows a unifrom distribution, in the limit $N_e arrow.r oo$ $c(omega) = 0$. For a normal distribution, on the average $c(omega) = N_e$ (this is the case of the electrons obtained from accelerators). Also, if the distribution has a very small variance (electrons are tightly confined in the bunch), we say that the distribution is #emph("microbunched"). In this case, $c(n dprod omega_u) = N_e^2$ where $n in NN^*$ is the number of transversal oscillations made by the bunch and $omega_u = (2 pi c) / lambda_u$ is the associated frequency. In reality, however, the microbunched scenario is hard to achieve.

For a real bunch, electrons do not have the same energy and momentum. The parameter that accounts for it is called #emph("emittance"), which is strongly linked to the volume occupied by the electrons in the 6D phase space $(vb(r)(t), vb(p)(t))$. A derived quantity is normalized emittance:

$ epsilon_a = sqrt(<Delta a^2> <Delta p^2 >- <Delta a Delta p_a>^2)/(m_e c), $

,where $Delta psi := psi - <psi>$ and $a in {x, y, z}$.

= Coumputational methods

In this chapter we are going to describe the main computaional tools in this study. We will start by introducing the PIC method and the particular implementation in FBPIC. Then we will discuss codes that are estimating the betatron radiation: the Synchrad code and the FBPIC implementaion that computes it on the fly. For the Byesian Optimization we will use Optimas.

== Particle-in-Cell codes

The Particle-in-Cell method (PIC) is the most widely used numerical approach to simulating laser-plasma phenomena because of its conceptual simplicity (the implementation is quite cumbersome) and its ability to capture the detailed spectra of accelerated particles. 

In this section, we will briefly introduce the method by following chapter 6.9 of @Gibbon2022. Essentially, these codes solve the self-consistent Vlasov-Maxwell equations (@Vlasov-eq + Maxwell). Solving this system numerically for real particles is not feasible, and because of this, we are reducing the dynamics of real particles to $N_("mp")$ #emph("macroparticles"). This approximation is valid because our system exhibits collective behavior. Each macro-particle can be thought of as a "brick" of the distribution function (see @vlasov-bricks):

$ f (vb(r), vb(p), t) := sum_i^(N_"mp") S_i (vb(r) - vb(r)_i (t)) delta (vb(p)- vb(p)_i (t)). $ <discrete-dist-function>

#figure(
  image("figures/vlasov-bricks.png", width:85%),
  caption: [Schematic representation of the #text(rgb("#00cc9a"))[distribution function] (a) and shattering of this function into numerical macroparticles (b). From @Pukhov2015],
) <vlasov-bricks>

$S_i (vb(r) - vb(r)_i (t))$ is the effective shape of the macroparticles, usually it is a Gaussian cenetered around $vb(r)_i$. By substituting @discrete-dist-function into @Vlasov-eq and first integrating over momentum and then for the position, the equations of motion for macroparticles are obtained:

$ dv(vb(p)_i, t) &= q_i/m_i (vb(E) + vb(v)_i cprod vb(B)) \
  dv(vb(r)_i, t) &= vb(p)_i / (gamma_i m_i), " " i=1, ..., N_("mp") $

Gathering the positions and velocities onto a grid, one can obtain the density and current needed to integrate Maxwell equations:

$ rho(vb(r)) &= sum_j q_j S (vb(r) - vb(r)_j) \
  vb(j)(vb(r)) &= sum_j q_j vb(v)_j S(vb(r) - vb(r)_j), " " j=1, ..., N_("cells"). $ <PIC-EM-sources>

Now, we can introduce the basic PIC loop (see @pic_step) as follows: starting from the initial positions and velocities of macroparticles, we can compute the EM sources with @PIC-EM-sources. With the computed sources, we integrate Maxwell equations to advance the EM fields. Now, we are interpolating new fields to the location of macroparticles and pushing them.

#figure(
  image("figures/pic_step.png", width: 58%),
  caption: [Schematic ilsutration of the PIC step. With #text(red)[red] there are quantities on the grid and with #text(blue)[blue] quantities on the position of macroparticles. Adapted from section 6.9 @Gibbon2022]
) <pic_step>