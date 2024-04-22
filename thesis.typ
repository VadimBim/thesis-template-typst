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

In the next sections, we will define plasma and how we can describe it. Then, we will limit our description to the physical scenario relevant to the LWFA case. After this, the dynamics of the accelerated bunch and its radiation will be briefly presented.

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

#grid(columns: (1fr, 1fr), math.equation(block: true, numbering: none, $ n_s (vb(r), t) = integral f_s dd(vb(p)^3) $), $ vb(u)_s (vb(r), t) = 1/n_s integral vb(v) f_s dd(vb(p)^3) $)

For a collisionless #footnote([If we want to include collisions in our model, we are inserting a #emph("collision integral") in the RHS of @Vlasov-eq]) plasma, the continuity equation in phase space is:

$ pdv(f_s, t) + vb(v) dprod grad_(vb(r)) f_s + q_s (vb(E) + vb(v) cprod vb(B)) dprod grad_(vb(p)) f_s = 0, $ <Vlasov-eq>

where $vb(v) = vb(p)/(m_s gamma_s) = (c vb(p))/sqrt(m_s^2 c^2 + vb(p)^2)$ is the velocity of particles in relativistic regime. The above equation is known as the Vlasov equation. Together with the Maxwell equations, it forms a #emph("self-consistent") system of equations that describes plasma dynamics. Where sources are given by:

#grid(columns: (1fr, 1fr), math.equation(block: true, numbering: none, $ rho = sum_s q_s n_s $), $ vb(j) = sum_s q_s n_s vb(u)_s $)

Where $q_s$ is the charge of species $s$. However, the Vlasov-Maxwell system is impossible to solve analytically except for a limited number of simple systems. In our case, we will employ a simpler model derived from kinetic theory (see section 2.2.2 @Macchi2013), the so-called #emph("fluid model") of plasma. this model assumes that each species with density $n_s$ and velocity $vb(u)_s$ behaves in a fluid-like manner. We will consider ion fluid stationary, and for the electrons we have:

$ &pdv(n_e, t) + div (n_e vb(u)_e) = 0 \
  &n_e m_e dv(vb(u)_e, t) = - n_e e (vb(E) + vb(u)_e cprod vb(B)) - grad P_e \
  &"Maxwell equations" \
  &dv(, t) (P_e n_e^(-gamma_e)) = 0 $ <fluid_model>

The last equation in @fluid_model is known as the closure condition, making our equations complete. $gamma_e = (2 + N)/N$ is the specific heat ratio, where $N$ is the number of degrees of freedom for the electron.

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

Taking this into account, @fluid_model_1D become:

$ pdv(n_1, t) + n_0 pdv(u_1, x) &= 0 \
  n_0 pdv(u_1, t) &= - e/m_e n_0 E_1 - 1/m_e pdv(P_1, x) \ 
  pdv(E_1, x) &= - e / epsilon_0 n_1 \
  P_1 &= 3 k_B T_e n_1 $ <fluid_model_1D_linearized>

For pressure, we assumed isothermal background electrons $P_0 = k_B T_e n_0$ (check @linearization_pressure for details). If we take time derivative of the first equation in @fluid_model_1D_linearized, write $pdv(u_1, t, s: \/)$ using second equation and eliminate $pdv(E_1, x, s: \/)$ and $P_1$ with third and fourth equation, we will get :

$ (pdv(, t, 2) - 3 v_("th")^2 pdv(, x, 2) + omega_p^2) n_1 = 0, $

where $omega_p$ is plasma frequency as defined in @plasma-lt-scale and $v_("th")^2 := k_B T_e slash m_e$ is the squared of electronic thermal velocity. If we look for plane wave solutions of the form $n_1 = n_0 exp(i(omega t - k x))$, differential operators become: $partial_t arrow.r i omega$ and $partial_x arrow.r -i k$, giving us in the end the Bohm-Gross dispersion relation:

$ omega^2 = omega_p^2 + 3 v^2_("th") k^2 $

== Betatron Radiation

The accelerated relativistic electrons obtained with LWFA can wiggle strongly. This wiggling causes electrons to radiate energy. The radiation emitted by an electron in the direction of observation $vb(n)$ (assuming that observation of the radiation is made #emph("far") from the electron) is given by (see eq. 14.65 from @Jackson1998-cw):

$ dd(cal(W), 2) / (dd(omega) dd(Omega)) = e^2/(4 pi^2 c) abs(integral_(-oo)^(+oo) dd(t) e^(i omega (t - (vb(n) dprod vb(r)(t))/c )) (vb(n) cprod [(vb(n) - vb(beta)) cprod dot(vb(beta))])/(1 - vb(beta) dprod vb(n))^2)^2. $

$dd(cal(W), 2) / (dd(omega) dd(Omega))$ is the energy radiated within a spectral band $dd(omega)$ centered on the frequency $omega$ in solid angle $dd(Omega)$ centered on the direction of observation $vb(n)$. $vb(beta)$ is the velocity of the electron normalized to the speed of light $c$, and $vb(r)(t)$ is the position of the electron function of time. Some important conlusions can be derived from the above equation:

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