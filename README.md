# Summer_project-df (2020/21 and 2021/22)
## Python GUI and Library for differential forms

This project aims at building a Python library for differential forms that gives an interactive visualization of the differential forms and their calculus such as exterior derivative, wedge product and Hodge operator, and helps understand their geometric meaning. The library will also be capable of plotting contravariant vector fields and visualizing different geometric structure such as gradient, divergence, curl and inner products between covariant and contravariant vectors. We aim at showing via graphical interface that differential forms are the natural object that should be used to describe fields and visualise their characteristics in different dimensions.

Some background: During 2020/21 academic year I (Moustafa) ran a Master's project that relied heavily on differential forms. My students struggled getting to grips with differential forms due to the lack of visual facility of differential forms analytically and numerically. Being familiar of a very old JAVA applet, [vector field analyser VFAII](https://math.la.asu.edu/~kawski/vfa2/vfa2sample.html), that gives excellent visualisation of vectors fields but with very limited 1-forms applications (stacks), I thought about "modernising" VFAII by writing it in Python with main focus on visualising and zooming on differential forms and doing their calculus.

Due to the ongoing pandemic, chances of academic travels are very low. So, I thought about running this idea as a Summer internship project as part of the School, of Physics and Astronomy at Nottingham University, Summer internship projects. That turns out to be a great plan as we managed to hire two excellent students for the project, Maciek and Sam.

The Summer 2021 project duration was 10 weeks from Monday 28/6/2021 to 5/9/2021. We met regularly, almost every day sometimes, and we had one long meeting (~ 2.5 hours) a week. The majority of meetings were online, but were very vibrant and productive; we bounced ideas, refined thoughts, optimised and fixed codes, discussed theory and derivations and planned for the next steps. Maciek and Sam did most of the coding in a fantastic fashion. 

This work gave birth to a Python GUI and a library, DFormpy. Both can do visualisation and calculus of vector fields on **R**<sup>2</sup> and visualisation and calculus of differential forms (0-form, 1-form and 2-form) on **R**<sup>2</sup>. The library has been published via a research paper. The paper gives detailed explanation of the library structure, plenty of examples and potential future implementations to the library. The publication can be found [here] (https://arxiv.org/abs/2201.10517)

In Summer 2022, we hired Vladimir as intern to work on my project. The Summer 2022 project aims at extending DFormpy functionality to deal with vectors and differential forms in **R**<sup>3</sup>.

