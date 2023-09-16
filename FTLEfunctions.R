# FTLE processing functions
# Updated 2023-07-05 to make different time label style, mid vs initial

library(animation)
library(patchwork)
library(raster)
library(PAMpal)
library(RCurl)
library(dplyr)
library(ggplot2)
library(scales)
library(PAMmisc)
library(ncdf4)


traceToDf <- function(nc, varname='ftle', backward=FALSE, days=2, timeStyle=c('initial', 'mid')) {
    timeStyle <- match.arg(timeStyle)
    trace <- stack(nc, varname=varname)
    allDf <- vector('list', length=length(trace@layers))
    for(i in seq_along(allDf)) {
        df <- as.data.frame(as(trace[[i]], 'SpatialPixelsDataFrame'))
        colnames(df) <- c('FTLE', 'Longitude', 'Latitude')
        df$name <- names(trace)[i]
        time <- as.POSIXct(names(trace)[i], format='X%Y.%m.%d.%H.%M.%S', tz='UTC')
        if(is.na(time)) {
            time <- as.POSIXct(names(trace)[i], format='X%Y.%m.%d', tz='UTC')
        }
        # forward we had half days, backward we subtract half days from initial
        if(timeStyle == 'mid') {
            time <- time + 24 * 3600 * days / 2 * ifelse(backward, -1, 1)
        }
        df$type <- ifelse(backward, 'Backward', 'Forward')
        df$time <- time
        allDf[[i]] <- df
    }
    allDf <- bind_rows(allDf)
    allDf$Longitude <- PAMmisc:::to180(allDf$Longitude)
    # allDf$time <- allDf$time + 24 * 3600 * days / 2
    allDf
}

plotFTLELoop <-function(trace, gps, gpsGroup=NULL,
                        xlim=NULL, ylim=NULL, zlim=NULL,
                        title=NULL, showPath=TRUE, progress=TRUE) {
    if(is.null(xlim)) {
        xlim <- range(c(gps$Longitude, trace$Longitude)) + c(-.01, .01)
    }
    if(is.null(ylim)) {
        ylim <- range(c(gps$Latitude, trace$Latitude)) + c(-.01, .01)
    }
    if(is.null(zlim)) {
        zlim <- range(trace$FTLE)
    }
    gps <- arrange(gps, UTC)
    if(!is.null(gpsGroup) &&
       !gpsGroup %in% colnames(gps)) {
        gpsGroup <- NULL
    }
    if(is.null(gpsGroup)) {
        starts <- gps[c(1, nrow(gps)), ]
        starts$location <- c('Start', 'End')
    } else {
        starts <- bind_rows(
            lapply(split(gps, gps[[gpsGroup]]), function(x) {
                x <- x[c(1, nrow(x)), ]
                x$location <- c('Start', 'End')
                x
            }
            )
        )
    }
    if(progress) {
        pb <- txtProgressBar(min=0, max=length(unique(trace$time)), style=3)
        ix <- 1
    }
    splitPlot <- 'type' %in% colnames(trace) && (length(unique(trace[['type']])) > 1)
    lapply(split(trace, trace$time), function(x) {
        if(splitPlot) {
            # browser()
            g <- plotFTLE(x[x$type == 'Forward', ], xlim, ylim, zlim,
                          title=NULL, gps, gpsGroup, starts=starts, showPath=showPath) +
                plotFTLE(x[x$type == 'Backward', ], xlim, ylim, zlim,
                         title=NULL, gps, gpsGroup, starts=starts, showPath=showPath) +
                plot_annotation(title)
        } else {
            g <- plotFTLE(x, xlim, ylim, zlim, title, gps, gpsGroup, starts, showPath=showPath)
        }
        print(g)
        if(progress) {
            setTxtProgressBar(pb, value=ix)
            ix <<- ix + 1
        }
    })
}

plotFTLEGIF <- function(trace, gps, gpsGroup='DriftName',
                        xlim=NULL, ylim=NULL, zlim=NULL,
                        title=NULL, progress=TRUE, height=800, width=1200, file, interval=.2,
                        showPath=TRUE) {
    saveGIF(
        plotFTLELoop(trace=trace, gps=gps, gpsGroup=gpsGroup,
                     xlim=xlim, ylim=ylim, zlim=zlim,
                     title=title, progress=progress, showPath=showPath),
        interval=interval, movie.name=file, ani.height=height, ani.width=width
    )
}

plotFTLE <- function(x, xlim=NULL, ylim=NULL, zlim=NULL, title=NULL, gps, gpsGroup=NULL, starts=NULL,
                     showPath=TRUE) {
    if(is.null(xlim)) {
        xlim <- range(c(gps$Longitude, x$Longitude)) + c(-.01, .01)
    }
    if(is.null(ylim)) {
        ylim <- range(c(gps$Latitude, x$Latitude)) + c(-.01, .01)
    }
    if(is.null(zlim)) {
        zlim <- range(x$FTLE)
    }
    if(!is.null(gpsGroup) &&
       !gpsGroup %in% colnames(gps)) {
        gpsGroup <- NULL
    }
    if(is.null(starts)) {
        if(is.null(gpsGroup)) {
            starts <- gps[c(1, nrow(gps)), ]
            starts$location <- c('Start', 'End')
        } else {
            starts <- bind_rows(
                lapply(split(gps, gps[[gpsGroup]]), function(s) {
                    s <- s[c(1, nrow(s)), ]
                    s$location <- c('Start', 'End')
                    s
                }
                )
            )
        }
    }
    if(is.null(gpsGroup)) {
        currentPoint <- data.frame(UTC = x$time[1])
        currentPoint <- addGps(currentPoint, gps, thresh=Inf)
    } else {
        currentPoint <- bind_rows(
            lapply(split(gps, gps[[gpsGroup]]), function(g) {
                cp <- data.frame(UTC=x$time[1])
                addGps(cp, g, thresh=Inf)
            })
        )
    }
    subtitle <- x$time[1]
    if('type' %in% colnames(x)) {
        uniType <- unique(x[['type']])
        if(length(uniType) == 1) {
            typeLabel <- switch(uniType,
                                'Forward' = 'Forward (Repelling) ',
                                'Backward' = 'Backward (Attracting) ')
            subtitle <- paste0(typeLabel, subtitle)
        }
    }
    g <- ggplot() +
        geom_tile(data=x, aes(x=Longitude, y=Latitude, fill=FTLE)) +
        xlim(xlim[1], xlim[2]) +
        ylim(ylim[1], ylim[2]) +
        # scale_fill_gradient(limits=zlim) +
        scale_fill_gradientn(
            # colors = viridis_pal(option = "plasma")(20),
            colors = colorRampPalette(c('white', 'seagreen', 'red', 'yellow'))(20),
            values = rescale(c(0,.5, 1, 2)),
            # for 24hr and 48hr integrations
            # values = rescale(c(0,.2,.25,.35, 0.4, 0.5,.6,.75, 1)),
            limits=c(0, 2),
            # for 96hr Integration
            # values = rescale(c(0,0.15,0.2,0.25,0.3,0.35, 0.4, 0.5)),
            # limits=c(0, .5),
            oob=scales::squish,
            na.value = "#FDE725FF") +
        ggtitle(title, subtitle=subtitle) +
        geom_point(data=starts, aes(x=Longitude, y=Latitude, col=location)) +
        geom_point(data=currentPoint, aes(x=Longitude, y=Latitude), col='black', size=3)
    if(showPath) {
        if(is.null(gpsGroup)) {
            g <- g + geom_path(data=gps, aes(x=Longitude, y=Latitude), col='black')
        } else {
            g <- g + geom_path(data=gps, aes(x=Longitude, y=Latitude, group=.data[[gpsGroup]]), col='black')
        }
    }
    g
}

gpsToHfradarDownload <- function(gps, days=2, buffer=1.5, name) {
    bbox <- c(max(gps$Latitude), min(gps$Latitude), min(gps$Longitude), max(gps$Longitude))
    bbox <- bbox + c(1, -1, -1, 1) * buffer
    # dayRange <- c(gps$UTC)), max(as.Date(gps$UTC))) + c(-1, 1) * days
    timeRange <- range(gps$UTC, na.rm=TRUE) + c(-1, 1) * days * 24 * 3600
    dayRange <- format(timeRange, format='%Y-%m-%d')
    hourRange <- format(timeRange, format='%H:%M:%S')
    hfnc <- hfRadar_Download(fname=name,
                             dStart = dayRange[1],
                             tStart=hourRange[1],
                             dEnd = dayRange[2],
                             tEnd=hourRange[2],
                             bbox = bbox)
    hfnc
}

transportFtpUpload <- function(file, pw, folder='ADRIFT/Raw') {
    ftpUp <- ftpUpload(file,
                       paste0('ftp://transport.me.berkeley.edu/',
                              folder,
                              '/',
                              basename(file)),
                       userpwd=pw)
    ftpUp
}

traceTimeSettings <- function(gps, days=2, increment=4) {
    dayRange <- c(min(as.Date(gps$UTC)), max(as.Date(gps$UTC))) + c(-1, 1) * days# * 24 * 3600
    dayDiff <- as.numeric(difftime(dayRange[2], dayRange[1], units='days'))
    (dayDiff - days) * (24 / increment)
}

traceSettings <-function(nc, days=2, increment=4) {
    nc <- nc_open(nc)
    on.exit(nc_close(nc))
    dim <- nc$dim
    names(dim) <- PAMmisc:::standardCoordNames(names(dim))
    nLon <- length(dim$Longitude$vals) * 10
    nLat <- length(dim$Latitude$vals) * 10
    times <- PAMmisc:::ncTimeToPosix(dim$UTC)
    hours <- as.numeric(difftime(max(times), min(times), units='hours')) - 24 * days
    nTimes <- floor(hours / increment)
    cat('Lat: ', nLat, ' Lon: ', nLon, ' Time: ', nTimes)
    invisible(c(nLat, nLon, nTimes))
}

