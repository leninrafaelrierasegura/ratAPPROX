get_sample <- function(n) {
  filename <- here::here("data_files", paste0("sample_n_", n, ".RData"))
  
  if (!file.exists(filename)) {
    set.seed(2)  # fixed seed → reproducible
    
    sample_data <- rnorm(n)  # <-- your sampling here
    
    save(sample_data, file = filename)
    cat("Generated and saved sample at:\n", filename, "\n")
  } else {
    cat("Loaded existing sample from:\n", filename, "\n")
  }
  
  load(filename)  # loads sample_data
  return(sample_data)
}