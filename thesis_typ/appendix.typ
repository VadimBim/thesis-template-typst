#counter(heading).update(1)
#set heading(numbering: "A.1")
#set math.equation(numbering: "1")

#import "@preview/physica:0.9.2": *
#import "@preview/cetz:0.2.2": canvas, plot


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

== Plasma's length-time scale <plasma-lt-scale>

In this section, we will quantatively define and describe two essential quantities in plasma physics: #emph("Debye length") $lambda_D$ and #emph("plasma frequency") $omega_p$.

Imagine we are inserting a positive charge $Q$ into a plasma:

#figure(
  image("../figures/big_charge.png",width:55%),
  caption: [The screening effect caused by a #text(rgb("#3333CC"))[positively] charged particle immersed into a plasma. Electrons are represented as #text(red)[-]])

In this case, electrons can screen the action of the "big" external charge. Distance at which this action can be felt is named Debye length $lambda_D$. We can derive it's value by considering the following conditions:

- Ions are immobile because $m_i >> m_e$
- Perturbation caused by the big charge is small: $(e Phi) / (k_B T) << 1$
- Thermal equilibrium is perserved: $p = n_e k_B T$

First, we will derive the density of electrons in this region. By using dynamical equilibrim:

$ vb(F) = -e vb(E) - 1/n_e grad p = vb(0), $

and thermal equilibrim condition we have:

$ &-e underbrace(vb(E), - grad Phi) = 1/n_e grad p = k_B T grad ln(n_e) arrow.r \
  &grad (k_B T ln(n_e) - e Phi) = vb(0) arrow.r \
  &n_e = n_0 exp((e Phi)/(k_B T)), $

where $n_0$ is the electron density before the charge was immersed and $Phi$ is the electric potential created by the charge. Recalling the Poisson's equation and using the assumption about small perturbations with quasi-neutrality ($n_i tilde.eq n_0$), we have:

$ laplacian Phi = - rho / epsilon_0 = - e/epsilon_0 (n_i - n_e) = - e/epsilon_0 n_0 underbracket((1 - exp((e Phi)/(k_B T))), tilde.eq (e Phi)/(k_B T)). $

Now using the spherical symmetry ($laplacian = 1/r^2  dv(, r) ( r^2 dv(Phi, r) )$) we arrive at the following second order differential equation:

$ dv(Phi, r, 2) + 2/r dv(Phi, r) - underbrace((n_0 e^2)/(epsilon_0 k_B T), := 1/lambda_D^2) Phi = 0, $ <Phi_ODE>

with the general solution: 

$ Phi(r) = c_1/r e^(-r/lambda_D) + c_2/r e^(#h(0.3em) r/lambda_D). $

By keeping only physically meaningful solution ($eval(Phi)_(r arrow.r oo) arrow.r 0$) and also knowing that we must return to classical Coulomb potential when there is no screening ($eval(Phi)_(lambda_D arrow.r oo) = Q/(4 pi epsilon_0) 1/r$) we arrive at:

$ Phi(r) = Q/(4 pi epsilon_0) e^(-r/lambda_D)/r . $ <screened_potential>

In @Phi_ODE we defined Debye length as follow:

$ lambda_D = sqrt((epsilon_0 k_B T)/(n_0 e^2)). $

We can observe that $lambda_D prop sqrt(T)$, which means that, from a phenomenological point of view, the higher the kinetic energy of the electrons is, the harder it is for them to screen the effect of the external charge. If the kinetic energy is small, electrons will be packed more densely around the external charge, thus increasing the screening effect ($lambda_D arrow.br$). Also, $lambda_D prop 1 / sqrt(n_0)$ states that the more electrons in the system, the better the screening of the external potential.

#figure(
canvas(length: 1.4cm, {
  plot.plot(size: (7, 5),
    x-tick-step: 1,
    y-tick-step: 1,
    x-max: 7,
    y-max: 6,
    x-min: 0,
    y-min: 0,
    x-label: [$r$ [a.u.]],
    y-label: [$Phi$ [a.u.]],
    legend: "legend.inner-north-east",
    {
      plot.add(
        domain: (0.01, 7),
        samples: 100,
        x => 1/x,
        label: [Coulomb potential]
        )
      plot.add(
        domain: (0.01, 7),
        samples: 100,
        x => 1/x * calc.exp(-x/2),
        label: [Screened potential $lambda_D = 2$]
        )
    })
}), caption: [Comparison between fast decaying screened potential (@screened_potential) and Coulomb potential. Bot $Phi$ and $r$ are in arbitrary units],)

Debye length and density of electrons determine if #emph("collective") effects will dominate over collisions between particles. The condition for this to happen is:

$ lambda_D >> "distance between particles", $

which is equivalent to the fact that the number of particles in the Debye sphere is big:

$ N_D = n_e (4 pi) / 3 lambda_D^3 >> 1. $

Now, we will define plasma frequency by looking at a simple model named #emph("capacitor model") (see @capacitor_model). Consider that we have a system of electrons and ions in equilibrium, and somehow, we display all electrons to the right by a small amount $x$. If we let the system evolve, electrons will oscillate (ions can be considered stationary because of the same assumption $m_i >> m_e$) with a specific frequency. To derive this frequency, we will start with the equation of motion for the electrons:

$ m_e dot.double(x) = -e E. $ <eom_capacitor>

The magnitude of the electric field will be given by the positive charge on the left surface $sigma$ in the approximation $sqrt("surface") >> x$ (see Example 2.2 from @Griffiths2023): 

$ E = sigma / epsilon_0 = (n_0 x)/epsilon_0. $ 

If we write $sigma$ in terms of volume charge density $sigma = n_0 x$, @eom_capacitor becomes:

$ dot.double(x) + omega^2_p x = 0, $

where we defined plasma frequency as:

$ omega_p = sqrt((n_0 e^2)/(m_e epsilon_0)). $

This quantity gives the characteristic time scale of the processes that occur in plasma.

#figure(
  image("../figures/capacitor_model.png", width: 85%),
  caption: [Capacitor model]
) <capacitor_model>

#pagebreak()

== Linearization of pressure equation <linearization_pressure>

#math.equation(block: true, numbering: none, $dv(, t) ((P_0 + P_1)/(n_0 + n_1)^3) &tilde.eq dv(, t) ((P_0 + P_1)/(n_0^3 + 3 n_0^2 n_1)) = 0 <=> \ 
  k_B T_e dv(, t) (1/(n_0^2 + 3 n_0 n_1)) &+ dv(, t) (P_1 / (n_0^3 + 3 n_0^2 n_1)) = 0 <=> \
  dv(, t) (P_1 / (n_0^3 + 3 n_0^2 n_1)) &= - k_B T_e dv(((n_0^2 + 3 n_0 n_1)^(-1)), (n_0^2 + 3 n_0 n_1)) dv((n_0^2 + 3 n_0 n_1), t) <=> \
  1/(n_0^3 + 3 n_0^2 n_1) dv(P_1, t) + P_1 dv(, t) (1/((n_0^3 + 3 n_0^2 n_1))) &= k_B T_e (3 n_0)/((n_0^2 + 3 n_0 n_1)^2) dv(n_1, t) <=> \ 
  1/(n_0^3 + 3 n_0^2 n_1) dv(P_1, t) - P_1 (3 n_0^2)/((n_0^3 + 3 n_0^2 n_1)^2) dv(n_1, t) &= (3 k_B T_e n_0)/((n_0^2 + 3 n_0 n_1)^2) dv(n_1, t) bar.v dprod (n_0^3 + 3 n_0^2 n_1)^2 => \
  (n_0^3 + 3 n_0^2 n_1) dv(P_1, t) - 3 n_0^2 P_1 dv(n_1, t) &= 3 n_0 k_B T_e dv(n_1, t) ((n_0^3 + 3 n_0^2 n_1)/(n_0^2 + 3 n_0 n_1))^2 <=> \
  dv(P_1, t) = 3 k_B T_e dv(n_1, t) <=> dv( ,t) (P_1 - 3 k_B T_e n_1) = 0 &<=> P_1 - 3 k_B T_e n_1 = "const." => \
  P_1 &= 3 k_B T_e n_1$)

In the transition from line 1 to 2 we used $P_0 = k_B T_e n_0$. For simplicity, we put constant equal to zero.