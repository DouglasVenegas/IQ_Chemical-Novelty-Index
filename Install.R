# install.R: Script para instalar paquetes necesarios en un entorno R

# Función para instalar y cargar paquetes de CRAN
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# Lista de paquetes de CRAN
cran_packages <- c(
  "repr", "IRdisplay", "evaluate", "crayon", "pbdZMQ", 
  "devtools", "uuid", "digest", "RColorBrewer", "tidyverse",
  "readxl", "rvest", "dplyr", "tidyr", "igraph", "visNetwork", 
  "readr", "ggplot2", "ggraph", "graphTweets", "writexl", 
  "httr", "proxy"
)

# Instalar y cargar paquetes de CRAN
install_and_load(cran_packages)

# Instalar IRkernel desde GitHub si no está instalado
if (!requireNamespace("IRkernel", quietly = TRUE)) {
  devtools::install_github("IRkernel/IRkernel")
  IRkernel::installspec(user = FALSE) # Registrar el kernel de R para Jupyter
}

# Mensaje de finalización
message("Todos los paquetes han sido instalados y están listos para usar.")
