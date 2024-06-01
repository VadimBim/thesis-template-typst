#let titlepage(
  title: "",
  degree: "",
  program: "",
  supervisor: "",
  advisors: (),
  author: "",
  startDate: none,
  submissionDate: none,
) = {
  set document(title: title, author: author)
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

  
  // --- Title Page ---
  v(1cm)
  align(center, image("../figures/antet_L.JPG", width: 100%))
align(center, image("../figures/eli_ifin_logo.jpeg", width: 85%))

v(10mm)

  align(center, text(font: sans-font, 1.0em, weight: 100, degree + "’s Thesis in " + program))

  v(3mm)
  
  line(length: 100%)
  smallcaps()[#align(center, text(font: sans-font, 1.3em, weight: 700, title))]
  line(length: 100%)


  v(2em)

  pad(
    top: 3em,
    right: 15%,
    left: 5%,
    grid(
      columns: 2,
      gutter: 1em,
      strong("Author: "), author,
    )
  )

  set align(right)

  pad(
    top: 5em,
    right: 1em,
    grid(
      columns: 2,
      gutter: 1em,
      strong("Supervisor: "), supervisor,
      strong("Faculty supervisor: "), advisors.join(", "),
    )
  )

  pagebreak()
}