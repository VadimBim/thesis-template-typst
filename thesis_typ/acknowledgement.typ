#let acknowledgement() = {
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

  set par(leading: 1em)

  
  // --- Acknowledgements ---
  align(left, text(font: sans-font, 2em, weight: 700,"Acknowledgements"))
  v(15mm)

  text[
    First of all, I want to thank my parents. I am deeply thankful for your unwavering support throughout my academic journey. 
    
    I would like to thank my advisor, Andrei Berceanu, for his insightful guidance, and invaluable expertise. I would also like to extend my heartfelt thanks to my faculty advisor, Mădălina Boca, for her continuous encouragement to pursue an academic career. Their perspectives have greatly enriched my work, and their support has been a cornerstone of my academic development.
    ]

}