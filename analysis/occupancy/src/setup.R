## add making data null for each model type
setup <- function(type, input) {
  if (type == "partner") {
    print("partner setup")
    input$monitors <-  c("phi.0",
                         "gam.0",
                         "phi.fra",
                         "gam.fra",
                         "phi.var",
                         "gam.var")
    vars <- c("X", "day", "day.2", "fra", "var.partner")
    input$data <- input$data[names(input$data)[grepl(paste(vars, collapse = "|"), names(input$data))]]
    remove.inits <- c("cnodf", "role")
    input$inits <- input$inits[!grepl(paste(remove.inits, collapse = "|"), names(input$inits))]
    names(input$inits)[grep("phi.var.partner", names(input$inits))] <- "phi.var"
    names(input$inits)[grep("gam.var.partner", names(input$inits))] <- "gam.var"
  }
  if (type == "all") {
    print("all setup")
    #input$data$degree <-NULL
    vars <- c("X", "day", "day.2", "fra", "cnodf", "role", "partner")
    remove.inits <- c("degree")
    input$inits <- input$inits[!grepl(paste(remove.inits, collapse = "|"), names(input$inits))]
  }
  if (type == "role" | type == "cnodf") {
    input$monitors <- c(
      "phi.0",
      "gam.0",
      "phi.fra",
      "gam.fra",
      "phi.var",
      "phi.mean",
      "gam.var",
      "gam.mean"
    )
    if (type == "role") {
      print("role setup")
      vars <- c("X", "day", "day.2", "fra", "role")
      input$data <- input$data[names(input$data)[grepl(paste(vars, collapse = "|"), names(input$data))]]
      remove.inits <- c("partner", "cnodf")
      input$inits <- input$inits[!grepl(paste(remove.inits, collapse = "|"), names(input$inits))]
      names(input$inits)[grep("phi.var.role", names(input$inits))] <- "phi.var"
      names(input$inits)[grep("phi.mean.role", names(input$inits))] <- "phi.mean"
      names(input$inits)[grep("gam.var.role", names(input$inits))] <- "gam.var"
      names(input$inits)[grep("gam.mean.role", names(input$inits))] <- "gam.mean"
    } else{
      print("cnodf setup")
      vars <- c("X", "day", "day.2", "fra", "cnodf")
      input$data <-
        input$data[names(input$data)[grepl(paste(vars, collapse = "|"), names(input$data))]]
      remove.inits <- c("partner", "role")
      input$inits <- input$inits[!grepl(paste(remove.inits, collapse = "|"), names(input$inits))]
      names(input$inits)[grep("phi.var.cnodf", names(input$inits))] <- "phi.var"
      names(input$inits)[grep("phi.mean.cnodf", names(input$inits))] <- "phi.mean"
      names(input$inits)[grep("gam.var.cnodf", names(input$inits))] <- "gam.var"
      names(input$inits)[grep("gam.mean.cnodf", names(input$inits))] <- "gam.mean"
    }
  }
  return(input)
}
