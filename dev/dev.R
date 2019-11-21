usethis::use_package("dplyr" ,min_version = "0.8.3")
usethis::use_package("tibble" ,min_version = "2.1.3")
usethis::use_package("taxize" ,min_version = "0.9.91")
usethis::use_package("taxizedb" ,min_version = "0.1.9.913")

devtools::document()


devtools::load_all()
