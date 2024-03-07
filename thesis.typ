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
    [Observe $f$ at $n$ initial points. $"Data" colon.eq cal(D) = {(bold(x)_1, y_1), ..., (bold(x)_n, y_n)}$],
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

More specifically, if $f$ is sampled from a GP, then a finite number function values $bold(f) colon.eq [f(bold(x)_1), f(bold(x)_2), ..., f(bold(x)_n)]^T$ follows a multivariate normal distribution:

$ bold(f) &tilde cal(N) (bold(mu), bold(K)) = (2 pi)^(-n slash 2) |bold(K)|^(-1 slash 2) exp ( - 1/2 (bold(f) - bold(mu))^T bold(K) (bold(f) - bold(mu))) $

, where $bold(mu)$ is the mean vector and $bold(K)$ is the $n times n$ covariance matrix given by the kernel $K_(i j) = k(bold(x)_i , bold(x)_j)$. Kernels have the property that points closer in input space are more strongly correlated: $ norm(bold(x) - bold(x)') < norm(bold(x) - bold(x)'') arrow.r.double k(bold(x), bold(x)') > k(bold(x), bold(x)'') $

If we want to make a new observation $f(bold(x)_(n+1)) colon.eq f_(n+1)$, by defintion of GP, it will come from the same probability distribution as $bold(f)$. Thus $P(f_(n+1) | bold(f))$ is obtained from marginalization of the underlying joint distribution $P([f_1, ... , f_n, f_(n+1)]^T)$, which is also a multivariate normal distribution. In conlusion. it can be shown @Rasmussen2005-ou that: // TODO put the derivation in the appendix.

$ &P(f_(n+1) | bold(f)) tilde cal(N) (mu_(n+1), sigma_(n+1)^(space 2)), "where" \
  &mu_(n+1) = k(x_(n+1), bold(x)) bold(K)^(-1) (bold(f) - bold(mu)) + mu(x_(n+1)) \
  &sigma_(n+1)^(space 2) = k(x_(n+1), x_(n+1)) - k(x_(n+1), bold(x)) bold(K)^(-1) k(bold(x), x_(n+1)) \
  &k(x_(n+1), bold(x)) colon.eq [k(x_(n+1), x_1), ..., k(x_(n+1), x_n)] = k(bold(x), x_(n+1))^T $

The covariance matrix should be symmetric and positive, which limits the number of possible functions that $k$ can take.