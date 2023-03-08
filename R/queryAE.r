get_request <-function(url){
  qr = URLencode(url)
  response <- GET( qr )
  
  if(status_code(response) != 200) {
    stop(
      paste(
        "Error running query. Received HTTP error code",
        status_code(response),
        "from server. Please try again later. If you continue to experience 
        problems please contact us 
        at https://www.ebi.ac.uk/about/contact/support/arrayexpress"
      )
    )
  }
  json_parsed <- fromJSON(txt = qr) 
  return(json_parsed)
}

queryAE = function(keywords = NULL, species = NULL){
  if(is.null(keywords) && is.null(species))
    stop("No keywords or species specified")
  
  baseURL <- "https://www.ebi.ac.uk/biostudies/api/v1/arrayexpress/search";
  page_size <- 100
  page=1
  
  query <- paste(baseURL,"?query=",keywords,"&organism=",species,"&pageSize=",
                 page_size,sep="")
  
  json_data <- get_request(paste(query,"&page=",page,sep=""))
  df_hits <- json_data$hits
  
  while (length(json_data$hits) > 0) {
    page <- page + 1
    json_data <- get_request(paste(query,"&page=",page,sep=""))
    if (length(json_data$hits) > 0) {
      df_hits <- rbind(df_hits, json_data$hits)
    }
  }
  return(df_hits)
}