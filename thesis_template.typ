#let project(
  title: "",
  degree: "",
  program: "",
  supervisor: "",
  advisors: (),
  author: "",
  startDate: none,
  submissionDate: none,
  body,
) = {
  set page(
    margin: (left: 30mm, right: 30mm, top: 40mm, bottom: 40mm),
    numbering: "1",
    number-align: center,
  )

  let body-font = "New Computer Modern"
  let sans-font = "New Computer Modern Sans"

  set text(
    font: body-font, 
    size: 12pt, 
    lang: "en"
  )
  
  show math.equation: set text(weight: 400)
  set math.equation(supplement: [Eq.])
  set math.equation(numbering: eqCounter => context {
    numbering(
      "(1.1.1.1)",
      ..counter(heading).get(),
      eqCounter
    )
  })

  // --- Headings ---
  show heading: h => {
    h
    // reset equation counter for each chapter
    counter(math.equation).update(0)
  }
  show heading: set block(below: 0.85em, above: 1.75em)
  show heading: set text(font: body-font)
  set heading(numbering: "1.1")
  // Reference first-level headings as "chapters"
  show ref: it => {
    let el = it.element
    if el != none and el.func() == heading and el.level == 1 {
      [Chapter ]
      numbering(
        el.numbering,
        ..counter(heading).at(el.location())
      )
    } else if el != none and el.func() == math.equation {
      el.supplement
      numbering(
        " (1.1.1.1)",
        ..counter(heading).at(el.location()),
        ..counter(math.equation).at(el.location())
      )
    } else {
      it
    }
  }

  // --- Paragraphs ---
  set par(leading: 1em)

  // --- Citations ---
  set cite(style: "alphanumeric")

  // --- Figures ---
  show figure: set text(size: 0.85em)
  set figure(supplement: [Fig.])
  
  // --- Table of Contents ---


  outline(
    title: {
      text(font: body-font, 1.5em, weight: 700, "Contents")
      v(15mm)
    },
    indent: 2em
  )
  
  show outline: set page(numbering: none)
  counter(page).update(1)

  v(2.4fr)
  pagebreak()


  // Main body.
  set par(justify: true, first-line-indent: 2em)

  body

  // Appendix.
  pagebreak()
  heading(numbering: none)[Appendix]
  include("thesis_typ/appendix.typ")

  pagebreak()
  bibliography("thesis.bib")
}
