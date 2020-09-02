# System requirements: see details in link https://docs.mosek.com/9.2/rmosek/install-interface.html
# 1. Windows
# Download Rtools (not an R package, link: https://cran.r-project.org/bin/windows/Rtools/). You will need to install the R toolset, the Cygwin DLLs, and the toolchain.
# Make sure you have mosek and the executables of Rtools on the PATH environment variable.
# 2. MacOS: Make sure you have Xcode installed.
# 3. Linux: On Ubuntu (and Debian) you may install the r-base-dev package.

installmosek <- function(targetDir = "~/bin"){
  if(!dir.exists(targetDir)){
    dir.create(targetDir)
  }
  if(!"Rmosek" %in% installed.packages()[,"Package"]) {
    install.packages("Rmosek")
  }
  switch(Sys.info()[['sysname']],
         Windows= {ostype <- "windows"},
         Linux  = {ostype <- "linux"},
         Darwin = {ostype <- "osx"})
  mosekUrl <- sprintf("https://download.mosek.com/stable/8.1.0.82/mosektools%s64x86.tar.bz2", ostype)
  download.file(mosekUrl, destfile = sprintf("%s/mosektools%s64x86.tar.bz2", targetDir, ostype))
  cmd <- sprintf("cd %s && tar -xvf mosektools%s64x86.tar.bz2", targetDir, ostype)
  system(cmd)
  Rmosek::mosek_attachbuilder(sprintf("%s/mosek/8/tools/platform/%s64x86/bin", targetDir, ostype))
  install.rmosek()
  if(!file.exists("~/mosek/mosek.lic")){
    message("Please visit: https://www.mosek.com/products/academic-licenses/ to get the licsense.
             Check email, put licsense file at: ~/mosek/mosek.lic")
  }
}
