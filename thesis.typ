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