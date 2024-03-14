#import "thesis_template.typ": *
#import "common/titlepage.typ": *
#import "thesis_typ/disclaimer.typ": *
#import "thesis_typ/acknowledgement.typ": *
#import "thesis_typ/abstract_en.typ": *
#import "common/metadata.typ": *

//Packages
#import "@preview/lovelace:0.2.0": *

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

- The input set is a hypter-rectangle $A = {bold(x) in RR^d | a_i lt.eq x_i lt.eq b_i}$ where $d$ is not large (typically, $d lt.eq 20$).
- The objective function is continuous.
- $f$ is "expensive to evaluate". This is especially true in the cases where objective is a result of a computationally demanding simulation, where each run can take hours.
- $f$ is a "black box" function, which means that we lack any knowledge about the structure of the function (linearity, concavity, periodicity, etc.). If we knew some of the properties, we could use more efficient algorithms that leverege this information.
- When we evaluate our objective, we observe only $f(bold(x))$ and no-first or second-derivatives.
- $f$ can be noisy, which is very relevant in the case where objective is obtained from an physical experiment.
- Also, BO is suited for cases where we want to find a global, rather than a local extremum.

The algorithm is composed of two main components: a statistical #emph("surrogate model") of the objective and an #emph("acqusition function"). The common choice for the surrogate is to use Gaussian Process (shortly GP) regression (@f_sim_gp) to model the function that we want to optimize. However, it can be anything that can make a prediction with the associated uncertainty about it.

$ f(bold(x)) tilde cal(G) cal(P) (mu(bold(x)), k(bold(x), bold(x)')) $ <f_sim_gp>

Informally, we can read @f_sim_gp as follow: the objective function follows a normal distribution over possible functions, determined by #emph("mean function") $mu(bold(x))$ and #emph("covariance function") $k(bold(x), bold(x)')$ (also known as the #emph("kernel")). 

The acqusition function $alpha(bold(x))$ is a heuristic that "tells" us where is the most promising point to evaluate our objective at the next iteration. The most promising point is given by the extremum of $alpha (bold(x))$:

$ bold(x)_n = "argmax" alpha_(n-1) (bold(x)) $

These functions are "cheaper" to optimize relative to the objective, so that traditional optimization techniques are employed. There are a lot of acqusition functions based on different assumptions. For example, a widely used acqusition function is Expected Improvemnt $upright(E)upright(I)(bold(x))$, as the name suggests this function returns the expected value to observe a better objective than observed so far. 

Basic Bayesian Optimization is a sequential algorithm, here the pseudo-code is shown (adapted from @Frazier2018-dq):

#show: setup-lovelace

#algorithm(
  caption: [Basic pseudo-code for Bayesian Optimization],
  pseudocode(
    [Place a Gaussian process prior on $f$],
    [Observe $f$ at $n$ initial points. $"Data" := cal(D) = {(bold(x)_1, y_1), ..., (bold(x)_n, y_n)}$],
    [*while* $n lt.eq N$ *do*], ind,
      [Update the posterior probability distribution on $f$ using $cal(D)$],
      [$bold(x)_(n + 1) = "argmax" alpha_(n)(bold(x))$],
      [Observe $y_(n+1) = f(bold(x)_(n+1))$],
      [Add $(bold(x)_(n+1), y_(n+1))$ to $cal(D)$], ded,
    [*end while*],
    [*return* best $y$]
  )
)

First $n$ points are usually generated using Sobol sequences @Sobol1967-nw to evenly sample from input space. Here $N$ is the number of available evaluations. 

In the next subchapters we disucss the components of BO in more details.

== Gaussian Processes <chapter-gp>

The formal defintion is as follow @Rasmussen2005-ou:

#rect(
  width: 100%,
  radius: 10%,
  stroke: 0.5pt,
  fill: rgb("#E6F9FF"),
)[#emph("A Gaussian process is a collection of random variables, any finite number of which have a joint Gaussian distribution")]

More specifically, if $f$ is sampled from a GP, then a finite number function values $bold(f) := [f(bold(x)_1), f(bold(x)_2), ..., f(bold(x)_n)]^T$ follows a multivariate normal distribution:

$ bold(f) &tilde cal(N) (bold(mu), bold(K)) = (2 pi)^(-n slash 2) |bold(K)|^(-1 slash 2) exp ( - 1/2 (bold(f) - bold(mu))^T bold(K)^(-1) (bold(f) - bold(mu))) $

, where $bold(mu)$ is the mean vector and $bold(K)$ is the $n times n$ covariance matrix given by the kernel $K_(i j) = k(bold(x)_i , bold(x)_j)$. Kernels have the property that points closer in input space are more strongly correlated: $ norm(bold(x) - bold(x)') < norm(bold(x) - bold(x)'') arrow.r.double k(bold(x), bold(x)') > k(bold(x), bold(x)'') $

If we want to make a new observation $f(bold(x)_(n+1)) := f_(n+1)$, by defintion of GP, it will come from the same probability distribution as $bold(f)$. Thus $P(f_(n+1) | bold(f))$ is obtained from marginalization (@marginalization) #footnote[All figures in @chapter-gp were generated using  http://www.infinitecuriosity.org/vizgp/.] of the underlying joint distribution $P([f_1, ... , f_n, f_(n+1)]^T)$, which is also a multivariate normal distribution. In conlusion. it can be shown @Rasmussen2005-ou that: // TODO put the derivation in the appendix.

$ &P(f_(n+1) | bold(f)) tilde cal(N) (mu_(n+1), sigma_(n+1)^(space 2)), "where" \
  &mu_(n+1) = k(x_(n+1), bold(x)) bold(K)^(-1) (bold(f) - bold(mu)) + mu(x_(n+1)) \
  &sigma_(n+1)^(space 2) = k(x_(n+1), x_(n+1)) - k(x_(n+1), bold(x)) bold(K)^(-1) k(bold(x), x_(n+1)) \
  &k(x_(n+1), bold(x)) := [k(x_(n+1), x_1), ..., k(x_(n+1), x_n)] = k(bold(x), x_(n+1))^T $ <posterior-dist>

#figure(
  image("figures/marginalization.png", width: 85%),
  caption: [The marginal distribution (#text(orange)[orange]) $P(f(2) | 3 "observations")$.],
) <marginalization>

In @marginalization we have an unkown 1D function, where we made 3 evaluations. Based on available data, we infered the probability distribution of $f(2)$.

The covariance matrix should be symmetric and positive, which limits the number of possible functions that $k$ can take. If we want to model smooth functions, the usual choice for the kernel is the Squared Exponential kernel (SE):

$ k(bold(x), bold(x')) &= a exp(- norm(bold(x) - bold(x'))^2), \ 
  norm(bold(x) - bold(x'))^2 &:= sum_(i=1)^d (x_i - x'_i) / (l_i^2) $ <SE_kernel>

$bold(theta) := [a, l_1, ..., l_d]$ are known as the #emph("hyperparameters") of the kernel, and they are responsible of the form of the sampled function. It is interesting to note that these hyperparameters are interpretable (see @samples_se). For example, if we take into consideration case $d=1$, we have only 2 hyperparameters $a$ and $l_1$. In this case, $a$ will influence the amplitude of the sampled functions, and $l_1$ the lengthscale (how fast our functions vary). 

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
  caption: [4 sampled functions from GPs with a SE kernel (bottom) and associated slices through the kernel $k(x_i, dot.c)$ (top) as a function of the second argument. On the left side $[a, l_1] = [1, 2]$ and on the right $[a, l_1^2] = [1.4, 0.08]$. ┈ are means, #box(rect(fill: rgb("#A2C6C6"), height: 8pt, width: 8pt, radius: 1pt), baseline: 5%) represents $plus.minus sigma$ and #box(rect(fill: rgb("#CADEDE"), height: 8pt, width: 8pt, radius: 1pt), baseline: 5%) $plus.minus 2 sigma$ confidence bands]
) <samples_se>

As we mentioned before, kernels are responsible for the structure of the modeled functions. Morevoer, we can combine different kernels via addition and multiplication to obtain functions with specific characteristics (see Chapter 2 of @Duvenaud2014).

#figure(
  image("figures/se_times_linear.png", width: 80%),
  caption: [Locally periodic 1D samples drawn from a GP with a  SE $times$ Per kernel. Where Per is a periodic kernel $:= a exp (-2 (sin^2(pi |x - x'| slash p ))/l^2)$]
)

For the mean function, the most common choice is a #emph("const.") value. It is possible to construct a mean that captures a specific trend of the function:

$ mu(bold(x)) = italic("const.") + sum_(i=1)^p beta_i psi_i (bold(x)), $

where $psi_i (bold(x))$ are known parametric functions, usually low-order polynomials. However, in the case of BO, we don't have paticular knowledge about the objective.

The modeled functions are very sensitive to the chosen hyperparameters, and the natural question is: How do we choose these hyperparameters such that sampled functions describe the best observed data? We will mention only Bayesian model selection, as the rest is out of scope of the current study (interested reader can consult Chapter 5 of @Rasmussen2005-ou for details). By Bayes rule, observed function values and hyperparameters are related as follow #footnote[Here we assumed that model is fixed by specifying form of the kernel and mean function, otherwise we should have written: $P(bold(theta), cal(M) | bold(f))$.]:

$ P(bold(theta) | bold(f)) = (P(bold(f) | bold(theta)) P(bold(theta))) / P(bold(f)) $

$P(bold(theta) | bold(f))$ is known as the #emph("posterior"), $P(bold(f) | bold(theta))$ is the #emph("likelihood") and $P(bold(f))$ is a normalization #emph("const.") named #emph("marginal likelihood") (a.k.a. evidence):

$ P(bold(f)) = integral P(bold(f) | bold(theta)) P(bold(theta)) d bold(theta) $

 Most of the times, this integral is not analytically tractable and analytical approximation or Monte carlo methods are used. As mentioned before, we want to maximize probability to see observed values given a set of hyperparameters :

$ hat(bold(theta)) = arg max P(bold(theta) | bold(f)) = arg max P(bold(f) | bold(theta)) P(bold(theta)), $ <MAP>

where we used the fact that $P(bold(f))$ doesn't depend on $bold(theta)$. In @MAP we are estimating hyperparameters by #emph("Maximum a posteriori") (MAP) estimate. An approximation can be obtained by assuming that $P(bold(theta))$ has a constant density over $bold(theta)$, such #emph("Maximum likelihood estimation") (MLE) (@MLE) is obtained. One can see that MAP is more computationally intensive, however it is a good choice when MLE gives unreasonable hyperparameters.

$ hat(bold(theta)) = arg max P(bold(f) | bold(theta)) $ <MLE>

Since GP assumes that function values are generated from a multivariate normal distribution, and taking the natural log we have:

$ ln P(bold(f) | bold(theta)) = -n/2 ln 2 pi - 1/2 ln |bold(K)| - 1/2 (bold(f) - bold(mu))^T bold(K)^(-1) (bold(f) - bold(mu)) $ <logP>

An important detail is that if our observations are noisy $bold(y) := bold(f) + bold(epsilon)$, we can model it by adding a diagonal matrix to the covariance matrix. In the case in which we assume that all obervations have the same variance (#emph("homoscedastic") noise), we simply add $sigma^2 I$ to $bold(K)$. When noise is different for each observation (#emph("heterorscedastic") noise), associated variances are added to the diagonal. Moreover, we can add correlations beetween noises to off-diagonal elements if known.

// TODO: derivation eq: 5.9 din @Rasmussen2005

Because $bold(theta)$ in @logP is nested inside the correlation matrix, deriving an analytical formula is not feasible. Therefore, we need to emply numerical optimization techinques to maximize $ln P(bold(f) | bold(theta))$.

== Acqusition functions <chapter-af>

In this chapter we will discuss in detail the most commonly used acqusition function, namely Expected Improvemnt. We will limit our analysis to 1D case and also will mention other functions used for #emph("standard") BO problems (check chapter 5 of @Frazier2018-dq for an introduction in Exotic BO).

Suppose that we performed a number of evaluations and $f(x^*)$ is the best value observed so far. Now we have one additional evaluation to perform, and we can perform it anywhere. We define improvemnt as $I (x) := max(f(x) - f(x^*), 0)$. While we would like to choose $x$ such that this improvement is large, $f(x)$ is unkown until after the evaluation. What we can do, however, is to take the expected value of this improvement and choose $x$ to maximize it. We define the expected improvement:

$ "EI" (x) := E [I (x)], $ <EI_definition>

where $"E"[dot.c]$ indicates the expectation taken under the posterior distribution . Expected improvement can be written in closed form (see @A_EI):

$ "EI" (x) &= (mu(x) - f(x^*)) (1 - Phi(z_0)) + sigma phi(z_0), \
  z_0 &= (f(x^*) - mu(x))/sigma(x) $ 

Another widely used acqusition function is #emph("Probability of Improvement"):

$ "PI" (x) := P(f(x) > f(x^*)). $ <PI_definition>

which can by simply evaluated with the same change of variables as in @A_EI, giving us :

$ "PI" (x) = Phi ((mu(x) - f(x^*)) / sigma(x)) $

The difference beetween PI and EI is that PI only considers the probability of improvement and not the expected magnitude of the new evaluation. The simplest acqusition function is #emph("Upper Confidence Bound"):

$ "UCB" := mu(x) + beta sigma(x), $

where $beta$ controls the exploration/exploitaion tradeoff. When $beta$ is small, the algorithm will explore regions with high $mu(x)$, on the contrary, when $beta$ is large, BO rewards the exploration of currently uncharted areas (high $sigma(x)$).