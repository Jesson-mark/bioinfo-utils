
check_files <- function(...){
    files <- list(...)
    # print(files)
    for(afile in files){
        # print(afile)
        if(! file.exists(afile)){
            stop('Error! File `', afile, '` does not exists! Please check it!')
        }
    }
}
