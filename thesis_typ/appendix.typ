#counter(heading).update(1)
#set heading(numbering: "A.1")
#set math.equation(numbering: "1")

#import "@preview/physica:0.9.2": *


== 1D Expected Improvemnt  <A_EI>

This derivation is adapted from #text("https://ekamperi.github.io/machine%20learning/2021/06/11/acquisition-functions.html", size: 0.55em).  

By definition we have:

$ "EI" (x) = "E" [I(x)] = integral_(-oo)^(+oo) I(x) cal(N) (mu(x), sigma^2(x)) d f(x) $

For simplicity, we droped index $n+1$ for $mu$ and $sigma$, but we must be aware that we are integrating over the posterior distribution given by previous observations. If we make the change of variables: $z = (f(x) - mu(x)) / sigma(x)$, we have:

$ &integral_(-oo)^(+oo) max(mu(x) + sigma(x) z - f(x^*), 0) phi(z) d z = \
&integral_(z_0)^(+oo) (mu(x) + sigma(x) z - f(x^*)) phi(z) d z = \
&(mu(x) - f(x^*)) integral_(z_0)^(+oo) phi(z) d z + sigma(x) integral_(z_0)^(+oo) z phi(z) d z = \
&(mu(x) - f(x^*))(1 - Phi(z_0)) + sigma(x) integral_(z_0)^(+oo) 1 / sqrt(2 pi) e^(-z^2 / 2) d (z^2/2) $ <EI_derivation>

, where $phi(z) = cal(N)(0, 1) = 1 / sqrt(2 pi) e^(-z^2 / 2)$ and $z_0 := (f(x^*) - mu(x)) / sigma(x)$. The second integral can be evaluated with a change of variables $u := z^2/2$. In the step from line 3 to 4 we used the definition of #emph("Cumulative Distribution Function") $:= Phi(z_0) = integral_(-oo)^(z_0) phi(z) d z $. The final result is: 

$ "EI"(x) = (mu(x) - f(x^*))(1 - Phi(z_0)) + sigma(x) phi(z_0) $ <EI_A>

We can see that EI is high when $mu(x) >> f(x^*)$, or when the uncertainty $sigma(x)$ is big. Another observation is that if we evaluate EI near an observed point ($sigma(x) tilde.eq 0$), then $"EI" tilde.eq 0 $ because $phi "and" (1 - Phi)$ at $+oo$ tends to zero.

If we want to control the exploration/exploitation trade-off, we can inject a parameter in @EI_A:

$ "EI" (x, xi) = (mu(x) - f(x^*) - xi )(1 - Phi(z'_0)) + sigma(x) phi(z'_0), $

where $z'_0 = (mu(x) - f(x^*) - xi )/sigma(x)$. In this way, we are pretending that the best current value is higher than it actually is #math.arrow.r more exploration.

== Plasma's lentgh-time scale <plasma-lt-scale>

In this section, we will quantatively define and describe two essential quantities in plasma physics: #emph("Debye length") $lambda_D$ and #emph("plasma frequency") $omega_p$.

Imagine we are inserting a positive charge into a plasma:

#figure(
  image("../figures/big_charge.png",width:55%),
  caption: [The screening effect caused by a #text(rgb("#3333CC"))[positively] charged particle immersed into a plasma. electrons are represented as #text(red)[-]])

In this case, surrounding charges can screen the action of the "big" external charge. Action at which this action can be felt is named Debye length. We can derive it's value by considering the following conditions:

- Ions are immobile because $m_i >> m_e$
- Perturbation caused by big charge is small: $(Delta U) / (k_B T) << 1$
- Thermal equilibrium is perserved: $p = n_e k_B T$

First, we will derive the density of electrons in this region. By using dynamical equilibrim:

$ vb(F) = -e vb(E) - 1/n_e grad p = vb(0), $

and thermal equilibrim condition we have:

$ &-e underbrace(vb(E), - grad Phi) = 1/n_e grad p = k_B T grad ln(n_e) arrow.r \
  &grad (k_B T ln(n_e) - e Phi) = vb(0) arrow.r \
  &n_e = n_0 exp((e Phi)/(k_B T)) $