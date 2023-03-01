library(Rcpp)


if( system.file(package="femR") != ""){
     system("rm src/*.o src/*.so")
     remove.packages("femR")
}

compileAttributes()

install.packages(".", type="source", repos=NULL)
