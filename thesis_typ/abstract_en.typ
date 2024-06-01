#let abstract_en() = {
  set page(
    margin: (left: 30mm, right: 30mm, top: 40mm, bottom: 40mm),
    numbering: none,
    number-align: center,
  )

  let body-font = "New Computer Modern"
  let sans-font = "New Computer Modern Sans"

  set text(
    font: body-font, 
    size: 12pt, 
    lang: "en"
  )

  set par(
    leading: 1em,
    justify: true
  )

  
  // --- Abstract ---
  v(1fr)
  align(center, text(font: body-font, 1em, weight: "semibold", "Abstract"))
  
  text[
    Laser–plasma interactions are of interest both for fundamental physics research (for example, reaching extreme electric fields to study QED phenomena) and as emerging technologies for industrial applications. In recent years, data-driven methods have proved to be efficient tools for exploring laser-plasma physics. 
    
    In this study, we applied a Bayesian machine-learning-based optimization algorithm to optimize the X-ray radiation obtained when a short pulse intense laser interacts with an underdense plasma. Computational simulations were carried out using Particle-in-Cell FBPIC code and the optimization library Optimas. Our studies showed that Bayesian Optimization could significantly reduce the number of iterations needed to find a promising physical scenario by tuning the relevant laser plasma input parameters. Future work can be done on incorporating physics knowledge into optimization scheme to reduce further the number of interations.
  ]
  
  v(1fr)
}