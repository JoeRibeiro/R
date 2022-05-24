rm(list=ls()) #clear environment, start afresh
MAINDIR<- "C:/Users/JR13/OneDrive - CEFAS/Fish_dataproduct_QSR_JR/"


#useful packages
years <- 1999:2011
quarters <- 1:4
library(icesDatras)

# Original getflexfile function in icesDATRAS library contains an error. it points to HH not ff
getFlexFile_fixed <- function(survey, year, quarter) {
  library(icesDatras)

  # check survey name
  if (!checkSurveyOK(survey)) return(FALSE)

  # check year
  if (!checkSurveyYearOK(survey, year, checksurvey = FALSE)) return(FALSE)

  # check quarter
  if (!checkSurveyYearQuarterOK(survey, year, quarter, checksurvey = FALSE, checkyear = FALSE)) return(FALSE)

  # read url and parse to data frame
  url <-
    sprintf(
      "https://datras.ices.dk/WebServices/DATRASWebService.asmx/getFlexFile?survey=%s&year=%i&quarter=%i",
      survey, year, quarter)
  print(url)
  out <- readDatras(url)
  out <- parseDatras(out)

  out
}

# Also this is not present in NAMESPACE
readDatras <- function(url) {
  # try downloading first:
  # create file name
  tmp <- tempfile()
  on.exit(unlink(tmp))

  # download file
  ret <-
    if (os.type("windows")) {
      download.file(url, destfile = tmp, quiet = TRUE)
    } else if (os.type("unix") & Sys.which("wget") != "") {
      download.file(url, destfile = tmp, quiet = TRUE, method = "wget")
    } else if (os.type("unix") & Sys.which("curl") != "") {
      download.file(url, destfile = tmp, quiet = TRUE, method = "curl")
    } else {
      127
    }

  # check return value
  if (ret == 0) {
    # scan lines
    scan(tmp, what = "", sep = "\n", quiet = TRUE)
  } else {
    message("Unable to download file so using slower method url().\n",
            "Try setting an appropriate value via\n\t",
            "options(download.file.method = ...)\n",
            "see ?download.file for more information.")
    # connect to url
    con <- url(url)
    on.exit(close(con))

    # scan lines
    scan(con, what = "", sep = "\n", quiet = TRUE)
  }
}



parseDatras <- function(x) {
  # parse using line and column separators
  type <- gsub(" *<ArrayOf(.*?) .*", "\\1", x[2])

  # convert any lazy teminated feilds to full feilds
  x <- gsub("^ *<(.*?) />$", "<\\1> NA </\\1>", x)
  starts <- grep(paste0("<", type, ">"), x)
  ends <- grep(paste0("</", type, ">"), x)
  ncol <- unique(ends[1] - starts[1]) - 1
  # drop everything we don't need
  x <- x[-c(1, 2, starts, ends, length(x))]

  # exit if no data is being returned
  if (length(x) == 0) return(NULL)

  # match content of first <tag>
  names_x <- gsub(" *<(.*?)>.*", "\\1", x[1:ncol])

  # delete all <tags>
  x <- gsub(" *<.*?>", "", x)
  # trim white space
  x <- trimws(x)

  # convert to data frame
  dim(x) <- c(ncol, length(x)/ncol)
  row.names(x) <- names_x
  x <- as.data.frame(t(x), stringsAsFactors = FALSE)

  # return data frame now if empty
  if (nrow(x) == 0) return(x)

  # DATRAS uses -9 and "" to indicate NA
  x[x == -9] <- NA
  x[x == ""] <- NA

  # simplify all columns except StatRec and AreaCode (so "45e6" does not become 45000000)
  x[!names(x) %in% c("StatRec", "AreaCode", "Ship")] <- simplify(x[!names(x) %in% c("StatRec", "AreaCode", "Ship")])

  x
}



# TODO - combine the check into readDatras - and do it at the download.file stage...
checkDatrasWebserviceOK <- function() {
  # return TRUE if web service is active, FALSE otherwise
  out <- readDatras("https://datras.ices.dk/WebServices/DATRASWebService.asmx/getSurveyList")

  # check server is not down by inspecting XML response for internal server error message
  if (grepl("Internal Server Error", out[1])) {
    warning("Web service failure: the server seems to be down, please try again later.")
    FALSE
  } else {
    TRUE
  }
}


simplify <- function(x) {
  # from Arni's toolbox
  # coerce object to simplest storage mode: factor > character > numeric > integer
  owarn <- options(warn = -1)
  on.exit(options(owarn))
  # list or data.frame
  if (is.list(x)) {
    for (i in seq_len(length(x)))
      x[[i]] <- simplify(x[[i]])
  }
  # matrix
  else if (is.matrix(x))
  {
    if (is.character(x) && sum(is.na(as.numeric(x))) == sum(is.na(x)))
      mode(x) <- "numeric"
    if (is.numeric(x))
    {
      y <- as.integer(x)
      if (sum(is.na(x)) == sum(is.na(y)) && all(x == y, na.rm = TRUE))
        mode(x) <- "integer"
    }
  }
  # vector
  else
  {
    if (is.factor(x))
      x <- as.character(x)
    if (is.character(x))
    {
      y <- as.numeric(x)
      if (sum(is.na(y)) == sum(is.na(x)))
        x <- y
    }
    if (is.numeric(x))
    {
      y <- as.integer(x)
      if (sum(is.na(x)) == sum(is.na(y)) && all(x == y, na.rm = TRUE))
        x <- y
    }
  }
  x
}


# returns TRUE if correct operating system is passed as an argument
os.type <- function (type = c("unix", "windows", "other"))
{
  type <- match.arg(type)
  if (type %in% c("windows", "unix")) {
    .Platform$OS.type == type
  } else {
    TRUE
  }
}

# Whilst not an ideal approach for downloading the data, it works. The library currently only supports one quarter/year at a time https://github.com/ices-tools-prod/icesDatras/issues/35
for(survey in getSurveyList()){
  i=0
  dl = NULL
  for(year in years){
    for(quarter in quarters){
        skip_to_next <- FALSE
        tryCatch(dl <- getFlexFile_fixed(survey,year,quarter), error = function(e) { skip_to_next <<- TRUE})
        if(skip_to_next | class(dl)!="data.frame") { next } else { 
          i=i+1 
          if(i==1){ #first file
            HH = dl} else { # append if not first successful download
              HH = rbind(HH,dl)
            }
        }
    }
  }
  if(i>0){write.csv(HH,paste0(MAINDIR,"HH_HL_download/HH/HH-",survey,".csv"))}
}

