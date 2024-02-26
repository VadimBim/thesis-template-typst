#let titlepage(
  title: "",
  titleGerman: "",
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
  align(center, image("../figures/antet_L.JPG", width: 120%))

  v(15mm)

  align(center, text(font: sans-font, 1.5em, weight: 100, degree + "’s Thesis in " + program))
  v(8mm)
  

  align(center, text(font: sans-font, 2em, weight: 700, title))

  pad(
    top: 3em,
    right: 15%,
    left: 15%,
    grid(
      columns: 2,
      gutter: 1em,
      strong("Author: "), author,
      strong("Supervisor: "), supervisor,
      strong("Advisors: "), advisors.join(", "),
    )
  )

  pagebreak()
}