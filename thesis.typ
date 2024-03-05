#import "thesis_template.typ": *
#import "common/titlepage.typ": *
#import "thesis_typ/disclaimer.typ": *
#import "thesis_typ/acknowledgement.typ": *
#import "thesis_typ/abstract_en.typ": *
#import "common/metadata.typ": *

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

= Bayesian Optimization

#emph("Bayesian Optimization") is a class of machine-learning-based algorithms that are designed to optimize expensive, black-box functions $f(bold(x))$ (from this point forward, we will refer to these functions as #emph("objective functions")). Usually, these algorithms are efficient when the objective has the follwing properties @Frazier2018-dq:

- The input set is a hypter-rectangle ${bold(x) in RR^d | a_i lt.eq x_i lt.eq b_i}$ where $d$ is not large (typicaly, $d lt.eq 20$).
- The objective functions is continuous.
- $f$ is "expensive to evaluate". This is especially true in the cases where objective is a result of a computationally demanding simulation, where each run can take hours.
- $f$ is a "black box" function, which means that we lack any knowledge about the structure of the function (linearity, concavity, periodicity, etc.). If we knew some of the properties, we could use more efficient algorithms that leverege this information.
- When we evaluate our obejective, we observe only $f(bold(x))$ and no-first or second-derivatives.
- $f$ can be noisy, which is very relevant in the case where objective is a obtained from an experiment.
- Also, Bayesian Optimization is suited for cases where we want to find a global, rather than a local extremum.

The algorithm is composed of two main components: a #emph("surrogate model") of the objective and an #emph("acqusition function"). For the surrogate, the common choice is to use Gaussian Process regression to model the function that we want to optimize.
However, it can be anything that can make a prediction with the associated uncertainty about it. The acqusition function $alpha(bold(x))$ is a heuristic that "tells" us where is the most promising point to evaluate our objective at. The most promising point is the extremum of $alpha (bold(x))$. These functions are "cheaper" to optimize relative to the objective, so that traditional optimization techniques are used. There are a lot of acqusition functions based on different assumptions. For example, a widely used acqusition function is Expected Improvemnt $E I (bold(x))$, as the name suggests this function returns the expected value to observe a better objective than observed so far. 

Basic Bayesian Optimization is a sequential algorithm, here the pseudo-code is shown:

In the next subchapters we disucss the components of Bayesian Optimization in more details.