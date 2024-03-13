#counter(heading).update(1)
#set heading(numbering: "A.1")
#set math.equation(numbering: "1")

== Expected Improvemnt  <A_EI>

$ "EI" (x) = "E" [I(x)] = integral_(-oo)^(+oo) I(x) cal(N) (mu(x), sigma^2(x)) d f(x) $

For simplicity, we droped index $n+1$ for $mu$ and $sigma$, but we must be aware that we are integrating over the posterior distribution given by previous observations. If we make the change of variables: $z = (f(x) - mu(x)) / sigma(x)$, we have:

$ &integral_(-oo)^(+oo) max(mu(x) + sigma(x) z - f(x^*), 0) phi(z) d z = \
&integral_(z_0)^(+oo) (mu(x) + sigma(x) z - f(x^*)) phi(z) d z = \
&(mu(x) - f(x^*)) integral_(z_0)^(+oo) phi(z) d z + sigma(x) integral_(z_0)^(+oo) z phi(z) d z = \
&(mu(x) - f(x^*))(1 - Phi(z_0)) + sigma(x) integral_(z_0)^(+oo) 1 / sqrt(2 pi) e^(-z^2 / 2) d (z^2/2) $ <EI_derivation>

, where $phi(z) = cal(N)(0, 1) = 1 / sqrt(2 pi) e^(-z^2 / 2)$ and $z_0 := (f(x^*) - mu(x)) / sigma(x)$. The second integral can be evaluated with a change of variables $u := z^2/2$. In the step from line 3 to 4 we used the definition of #emph("Cumulative Distribution Function") $:= Phi(z_0) = integral_(-oo)^(z_0) phi(z) d z $. The final result is: 

$ "EI"(x) = (mu(x) - f(x^*))(1 - Phi(z_0)) + sigma(x) phi(z_0) $ <EI_A>

We can see that EI is high when $mu(x) >> f(x^*)$, or when the uncertainty $sigma(x)$ is big. Another observation is that if we evaluate EI near an observed point ($sigma(x) tilde.eq 0$), then $"EI" tilde.eq 0 $ because $phi "and" (1 - Phi)$ at $+oo$ tends to zero.

If we want to contral the exploration/explotation trade-off, we can inject a parameter in @EI_A:

$ "EI" (x, xi) = (mu(x) - f(x^*) - xi )(1 - Phi(z'_0)) + sigma(x) phi(z'_0), $

where $z'_0 = (mu(x) - f(x^*) - xi )/sigma(x)$. In this way, we are pretending that the best current values is higher that it actually is #math.arrow.r more exploration.