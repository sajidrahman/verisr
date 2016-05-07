#' get matrix for CSV
#'
#' given a verisr object, this will create a matrix of the data
#' suitable for saving off as a CSV.
#'
#' @param veris a verisr object
#' @export
veris2csv <- function(veris) {
  vnames <- getvnames(veris)
  # setup venum with field names to represent as categories
  # the names must be passable to getenumlist
  enumfields <- c("^actor", "^action", "availability.variety", "integrity.variety")
  venum <- vnames[grep(paste(enumfields, collapse='|'), vnames)]
  venum <- venum[grep("cve|name|notes|country|industry", venum, invert=T)]
  venum <- c(venum, "actor", "action")
  # vtext are fields to just dump as text
  vtext <- vnames[grep(paste(venum, collapse='|'), vnames, invert=T)]
  # ignore the one offs for now (do these below)
  oneoff <- c("attribute.confidentiality.data\\.",
              "asset.assets",
              "impact.loss",
              "schema_version")
  # pull them out of vtext
  vtext <- vtext[grep(paste(oneoff, collapse='|'), vtext, invert=T)]
  # pull the text into a results for text field
  rtext <- sapply(vtext, function(x) {
    ## victim will be messed up
    rval <- sapply(getenumlist(veris, x), function(vrow) { unlist(vrow)[1] })
    rval[is.na(rval)] <- ""
    rval
  })
  # now do the categorical values
  renum <- do.call(cbind, sapply(venum, function(x) {
    # create a matrix for this one venum name, one column for each value pulled
    raw <- getenumlist(veris, x)
    rawname <- sort(unique(unlist(raw)))
    outmat <- matrix(data=c(0), nrow=length(veris), ncol=length(rawname))
    colnames(outmat) <- paste(x, rawname, sep='.')
    for(i in seq(raw)) {
      outmat[i, which(rawname %in% raw[[i]])] <- 1
    }
    outmat
  }))
  ## unique asset categories
  my.assets <- c("asset.assets", "asset.assets.variety")
  assets <- do.call(cbind, lapply(my.assets, function(x) {
    raw <- getenumlist(veris, x)
    rawname <- sort(unique(unlist(raw)))
    asset <- matrix(data=c(0), nrow=length(veris), ncol=length(rawname))
    colnames(asset) <- paste(x, rawname, sep='.')
    for(i in seq(raw)) {
      thisval <- table(unlist(raw[[i]]))
      asset[i, which(rawname %in% names(thisval))] <- as.vector(thisval)
    }
    asset
  }))
  #sortvmatrix(cbind(rtext, renum, assets))
  sortvmatrix(cbind(data.frame(rtext), data.frame(renum), data.frame(assets)))
}

#' Internal: sort veris columns
#'
#' This will sort the columns of a matrix so Wade is happy.
#'
#' @param x a matrix where columns are veris fields
sortvmatrix <- function(x) {
  sortorder <- c("incident_id", "source_id", "reference", "security_incident",
                 "confidence", "summary", "related_incident",
                 "notes", "victim", "actor.ext", "actor.int", "actor.par",
                 "actor.unk",
                 "action.mal", "action.hack", "action.soc", "action.mis",
                 "action.phy", "action.err", "action.env", "action.unk",
                 "asset",
                 "attribute.conf", "attribute.int", "attribute.avail",
                 "timeline", "discovery", "target", "control", "correct",
                 "cost", "ioc", "impact", "plus")
  unsorted <- colnames(x)
  sortnames <- NULL
  for(i in seq_along(sortorder)) {
    sortnames <- c(sortnames, unsorted[grep(paste("^", sortorder[i], sep=''), unsorted, perl=T)])
  }
  missed <- grep(paste("^", sortnames, sep="", collapse="|"), unsorted, perl=T, invert=T)
  if(length(missed)) {
    sortnames <- c(sortnames, unsorted[missed])
  }
  x[ , sortnames]
}

#' convert veris object to a matrix
#'
#' given a verisr object, this will create a matrix of the data
#' suitable for further PCA analysis
#'
#' @param veris a verisr object
#' @param industry either 2 or 3 digits for industry
#' @param unknown logical, whether to include Unknown/Other fields
#' @export
veris2matrix <- function(veris, industry=2, unknown=F) {
  vnames <- getvnames(veris)
  enumfields <- c("actor", "action")
  venum <- vnames[grep(paste('^',enumfields, sep='', collapse='|'), vnames)]
  venum <- venum[grep("cve|name|notes|country|industry", venum, invert=T)]
  venum <- c(venum, vnames[grep('timeline.*unit', vnames)] )
  if (industry==2) {
    venum <- c(venum, "victim.industry2")
  } else {
    venum <- c(venum, "victim.industry3")
  }
  venum <- c(venum, "actor", "action",
             "victim.employee_count",
             "security_incident", "asset.assets", "asset.assets.variety", "asset.cloud",
             "asset.hosting", "asset.management", "asset.ownership",
             "attribute.confidentiality.data.variety", "attribute.confidentiality.data_disclosure",
             "discovery_method", "targeted",
             "attribute.integrity.variety", "attribute.availability.variety")
  # skipping "discovery_method", "targeted"
  # fixed that, not skipping discovery_method or targeted
  # now do the categorical values
  renum <- do.call(cbind, sapply(venum, function(x) {
    # create a matrix for this one venum name, one column for each value pulled
    raw <- getenumlist(veris, x)
    rawname <- sort(unique(unlist(raw)))
    if (length(rawname)==0) {
      return(NULL)
    }
    outmat <- matrix(data=c(0), nrow=length(veris), ncol=length(rawname))
    colnames(outmat) <- paste(x, rawname, sep='.')
    for(i in seq(raw)) {
      outmat[i, which(rawname %in% raw[[i]])] <- 1
    }
    outmat
  }))
  final.matrix <- sortvmatrix(renum)
  # strip out unknowns, it's the right thing to do.
  # careful here, also stripping out "NA" which may not be appropriate
  if(!unknown) {
    final.matrix <- final.matrix[ ,grep("Unknown|unknown|00|NA", colnames(final.matrix), invert=T)]
  }
  final.matrix[ ,colSums(final.matrix)!=nrow(final.matrix)]
}

#' compact rows for analysis
#'
#' given a verisr object and a list of features, this will create a matrix of the data
#' suitable for further PCA, hclust or other dist() analysis.
#'
#' @param pca a pca matrix object
#' @param feature a vector of features to fold on
#' @param min if min is > 1, this will check for a minimum amount of records for each feature,
#' if it doesn't have the min, it will be set to zeros and cleaned if clean is set to TRUE
#' @param clean if TRUE, this will look for anything that will mess up a PCA or MDS function and remove it.
#' @export
foldmatrix <- function(pca, feature, min=1, clean=FALSE) {
  folded <- do.call(rbind, lapply(feature, function(x) {
    # pull slice of rows where that incident has feature x
    incidents <- which(pca[ ,x]==1)
    if (length(incidents) < max(2, min)) {
      ret <- pca[1, ]
      ret[ret>0] <- 0
    } else {
      pcas <- pca[incidents, ]
      ret <- colMeans(pcas)
    }
    ret
  }))
  rownames(folded) <- feature
  if (clean) {
    # removing zero-sum rows (addressing "min" above)
    allsums <- rowSums(folded)
    # and columns need non-zero variance
    allvars <- apply(folded, 2, var)
    folded <- folded[which(allsums>0), which(allvars>0)]
  }
  folded
}

#' fold rows of matrix for further analysis
#'
#' given a verisr object and a list of features, this will create a matrix of the data
#' suitable for further PCA analysis.
#'
#' @param veris a verisr object
#' @param filter an optional filter to apply
#' @export
pcamatrix <- function(veris, filter=NULL) {
  if (!is.null(filter)) {
    if (is.logical(filter)) {
      veris <- veris[filter]
    } else {
      stop("filter is not logical vector")
    }
  }
  # get the raw PCA matrix
  rawmatrix <- veris2matrix(veris)
  allnames <- colSums(rawmatrix)
  cat("Folding", nrow(rawmatrix), "rows by", ncol(rawmatrix), "cols.\n")
  # have at least 1% or 5, whichever is greater
  venum <- names(allnames[which(allnames>(max(nrow(rawmatrix)*0.01, 8)))])

  feature <- venum[grep("^actor|^action|^asset.assets.variety|data.variety", venum)]

  pca.ready <- foldmatrix(rawmatrix, feature)
  # want foldvmatrix to take in a list of features or a list of filters
  # for now will take on just the feature list.

  pca.ready
}

#' count incidents in clusters
#'
#' given a verisr object and a list of labels, this will create a data.frame of the data
#' suitable for plotting.
#'
#' @param veris a verisr object
#' @param label a list of labels to count
#' @param filter an optional filter to apply
#' @export
countMatrix <- function(veris, label, filter=NULL) {
  if (!is.null(filter)) {
    if (is.logical(filter)) {
      veris <- veris[filter]
    } else {
      stop("filter is not logical vector")
    }
  }
  rawmatrix <- veris2matrix(veris)

  folded <- do.call(rbind, lapply(seq_along(label), function(cluster) {
    do.call(rbind, lapply(label[[cluster]], function(x) {
      data.frame(cluster=cluster, enum=x, count=sum(rawmatrix[ ,x]==1))
    }))
  }))
  folded
}

#'following functions copied from csv.R file'
#' Internal: Get fields names
#'
#' This will grab all the field names from a veris object
#'
#' @param veris a verisr object
getvnames <- function(veris) {
  sort(unique(unlist(getvnamelong(veris))))
}

#' Internal: Get fields names using a long method
#'
#' This will grab all the field names from a veris object.
#'
#' @param veris a verisr object
#' @param curname used to maintain state internally
getvnamelong <- function(veris, curname = NULL) {
  if(is.null(names(veris)) & length(veris)) {
    return(lapply(veris, getvnamelong, curname))
  }
  realname <- function(n) {
    if (is.null(curname)) n else paste0(curname, ".", n)
  }
  allthenames <- lapply(names(veris), function(x) {
    if (mode(veris[[x]]) %in% c("character", "numeric", "logical")) {
      ret <- realname(x)
    } else {
      ret <- getvnamelong(veris[[x]], realname(x))
    }
    ret
  })
  allthenames
}

#' Internal: sort veris columns
#'
#' This will sort the columns of a matrix so Wade is happy.
#'
#' @param x a matrix where columns are veris fields
#sortvmatrix <- function(x) {
#  sortorder <- c("incident_id", "source_id", "reference", "security_incident",
#                 "confidence", "summary", "related_incident",
#                 "notes", "victim", "actor.ext", "actor.int", "actor.par",
#                 "actor.unk",
#                 "action.mal", "action.hack", "action.soc", "action.mis",
#                 "action.phy", "action.err", "action.env", "action.unk",
#                 "asset",
#                 "attribute.conf", "attribute.int", "attribute.avail",
#                 "timeline", "discovery", "target", "control", "correct",
#                 "cost", "ioc", "impact", "plus")
#  unsorted <- colnames(x)
#  sortnames <- NULL
#  for(i in seq_along(sortorder)) {
#    sortnames <- c(sortnames, unsorted[grep(paste("^", sortorder[i], sep=''), unsorted, perl=T)])
#  }
#  missed <- grep(paste("^", sortnames, sep="", collapse="|"), unsorted, perl=T, invert=T)
#  if(length(missed)) {
#    sortnames <- c(sortnames, unsorted[missed])
#  }
#  x[ , sortnames]
#}

#' return a list matching vcdb ordering and length with requested object
#'
#' This will iterate through the veris object and return
#' a list of matches.  This is intented to maintain the orginal
#' indexes of the veris object so further manipulation can be done.
#'
#' Note: Can do a special "industryN" request and it will chop
#' off the industry at the N value or return same length of zeros
#' if it isn't long enough.
#'
#' @param veris a verisr object
#' @param enum the field to count
getenumlist <- function(veris, enum) {
  # if the veris object has null names and yet length
  # it is an array, and we simply want to step into
  # and through it.  The top level veris object
  # is an array, as is things like victim and assets
  # and data variety
  if(is.null(names(veris)) & length(veris)) {
    return(lapply(veris, getenumlist, enum))
  }
  # now we are in the meat of the function
  # and we should have either a full slice
  # or a partial slice of a veris incident

  # look at the enum passed in, want to
  # grab the first ("tag") and concatenate the rest
  vars <- unlist(strsplit(enum, "[.]"))
  tag <- vars[1]
  therest <- paste(vars[-1], collapse='.')
  # if the veris object is null at "tag", return NA
  if (is.null(veris[[tag]])) {
    retval <- NA
  } else if (therest == "") {
    # else if we are at the end of our enum, return the value?
    if (length(veris[[tag]])==0) {
      retval <- NA
    }
    # if we have names return those
    # it's an easy way to count actions, actors, etc.
    these.names <- names(veris[[tag]])

    #cat("the rest is blank, names:", these.names, "null:", is.null(these.names), "\n")
    if (tag=="assets") {
      assetmap <- c("S"="Server", "N"="Network", "U"="User Dev", "M"="Media",
                    "P"="Person", "T"="Kiosk/Term", "Unknown"="Unknown")
      retval <- unique(unlist(sapply(veris[[tag]], function(asset) {
        myasset <- ifelse(asset$variety=="Unknown", "Unknown", substr(asset$variety, 1, 1))
        myamount <- 1 # not counting more than one here
        # myamount <- ifelse(is.null(asset$amount), 1, asset$amount)
        rep(assetmap[[myasset]], myamount)

      })))
    } else if (!is.null(these.names)) {
      retval <- these.names
      # note to self, this is causing the getMatrix functions
      # to return NA, as the names of the return vector are blank
      #    } else if (is.null(these.names)) {
      #      retval <- NA
    } else {
      retval <- veris[[tag]]
    }
  } else {
    # else we need to continue to "drill down" into the veris object
    # with the rest of the enum being quieried
    # passing it back to self so it can parse through arrays
    # and use the same logic to continue parsing
    #
    # but before we do, let's check for some unique variables
    # like "industry*" where * is a length to chop
    if (grepl("^industry\\d$", therest, perl=T)) {
      # figure out the length of industry to return
      ind.len <- substr(therest, 9, 9)
      retval <- getenumlist(veris[[tag]], "industry")
      retval <- lapply(retval, function(x) {
        i <- substr(x, 1, ind.len)
        ifelse(nchar(i)==ind.len, i, paste(rep("0", ind.len), collapse=""))
      })
    } else if (tag=="assets" & therest=="variety") {
      retval <- unlist(sapply(veris[[tag]], getVarietyAmount))
    } else if (tag=="data" & therest=="variety") {
      retval <- unlist(getenumlist(veris[[tag]], therest))
    } else if (tag=="loss" & therest=="variety") {
      # TODO: impact could use attention
      retval <- unlist(sapply(veris[[tag]], getVarietyAmount))
    } else {
      retval <- getenumlist(veris[[tag]], therest)
    }
  }
  retval
}

#' Internal: Expand variety
#'
#' This will expand all of the "variety" fields by the amount
#' specified in the "amount" field of the same level object.
#'
#' @param x a slice of a veris a verisr object
getVarietyAmount <- function(x) {
  variety <- x[['variety']]
  #if ('amount' %in% names(x)) {
  #  amount <- ifelse(x[['amount']]>1, x[['amount']], 1)
  #  variety <- rep(variety, amount)
  #}
  variety
}
