#!/usr/bin/env Rscript

## Usage Script:

## Rscript circos_numts.R <directory containing final NUMTs> <ref_headers.txt> <reference fasta index file> <MT Length>

required <- c("tidyverse", "RColorBrewer", "circlize", "svglite", "optparse", "stringr")

sapply(required, function(x) {
    if(!require(x, character.only = TRUE)) {
        install.packages(x)
        suppressPackageStartupMessages(library(x, character.only=TRUE))
    } else {
        suppressPackageStartupMessages(library(x, character.only=TRUE))
    }
})

pdf(NULL)

option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="Input File", metavar="character"),
    make_option(c("-s", "--svg"), type="character", default=NULL,
                help="Output SVG file name", metavar="character"),
    make_option(c("-p", "--png"), type="character", default=NULL,
                help="Output PNG file name", metavar="character"),
    make_option(c("-r", "--ref"), type="character", default=NULL,
                help="Nuclear Reference Headers List", metavar="character"),
    make_option(c("-f", "--fai"), type="character", default=NULL,
                help="Nuclear Reference Fasta Index", metavar="character"),
    make_option(c("-l", "--length"), type="numeric", default=NULL,
                help="Mitochondrial Reference Length", metavar="numeric")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 6 ) {
    print_help(opt_parser)
    stop("Please provide all the required inputs.", call.=FALSE)
}

input <- opt$input

out_png <- opt$png

out_svg <- opt$svg

required_ref_headers <- opt$ref

ref_index <- opt$fai

mt_length <- opt$length

mt_length <- as.numeric(mt_length)

data <- read.delim(file = input, header=FALSE)

colnames(data) <- c("Chromosome", "Start", "MT_Header", "MT_Start", "MT_End", "Length")

data$End <- data$Start+1

data <- data %>% dplyr::select(c(Chromosome, Start, End, MT_Header, MT_Start, MT_End, Length))

data <- data %>% mutate(Type=ifelse(abs(((MT_End-MT_Start)-Length))>abs(((mt_length+MT_Start)-MT_End)-Length), "Control", "Non-Control"))

data_control <- data %>%
                    filter(Type=="Control") %>%
                    dplyr::select(c(Chromosome, Start, End, MT_Header, MT_Start, MT_End)) %>%
                    pivot_longer(cols=-c(Chromosome, Start, End, MT_Header), names_to="Type", values_to="MT_Coordinates") %>%
                    mutate(MT_Start = ifelse(Type=="MT_Start", mt_length+1, MT_Coordinates), MT_End = ifelse(Type=="MT_Start", MT_Coordinates+mt_length, mt_length)) %>%
                    dplyr::select(c(Chromosome, Start, End, MT_Header, MT_Start, MT_End))

data_non_control <- data %>%
                        filter(Type=="Non-Control") %>%
                        dplyr::select(c(Chromosome, Start, End, MT_Header, MT_Start, MT_End))

final_data <- rbind(data_control, data_non_control)

final_data <- final_data %>% mutate(MT_Start_1 = ifelse(MT_Start>MT_End, MT_End, MT_Start),
                        MT_End_1 = ifelse(MT_End > MT_Start, MT_End, MT_Start)) %>%
                        dplyr::select(c(Chromosome, Start, End, MT_Header, MT_Start_1, MT_End_1))

colnames(final_data)[5:6] <- c("MT_Start", "MT_End")

final_data <- final_data %>% mutate(MT_Header = ifelse(MT_Start>mt_length & MT_End>mt_length, "Concat_2", "Concat_1"))

final_data <- as.data.frame(lapply(final_data, unlist))

final_data$Chromosome <- gsub("[C|c]hr", "", final_data$Chromosome)

index_data <- read.delim(file = ref_index, header = FALSE)

index_data <- index_data[,1:2]

colnames(index_data) <- c("Chromosome", "Max")

index_data$Chromosome <- as.character(index_data$Chromosome)

target_headers <- read.delim(file = required_ref_headers, header = FALSE)

colnames(target_headers) <- "Chromosome"

target_headers$Chromosome <- stringr::str_sort(as.character(target_headers$Chromosome), numeric=TRUE)

chromosome_data <- inner_join(target_headers, index_data, by = "Chromosome") %>% dplyr::select(c("Chromosome", "Max"))

chromosome_data$Chromosome <- gsub("[C|c]hr", "", chromosome_data$Chromosome)

chromosome_data$Min <- rep(1,nrow(chromosome_data))

chromosome_data <- chromosome_data %>% dplyr::select(c("Chromosome", "Min", "Max"))

sector <- c(chromosome_data$Chromosome, "Concat_1", "Concat_2")

sector_xlim <- chromosome_data[,2:3]

sector_xlim <- rbind(sector_xlim, c(1, mt_length), c(mt_length+1, 2*mt_length))

mt_col <- c("#5ab4ac", "#C5C5C5")

names(mt_col) <- c("Concat_1", "Concat_2")

if(nrow(chromosome_data)<=12) {
    nucl_col <- brewer.pal(nrow(chromosome_data), "Set3")
    names(nucl_col) <- as.character(chromosome_data$Chromosome)
} else if(nrow(chromosome_data)>12 && nrow(chromosome_data)<=20) {
    nucl_col <- c(brewer.pal(12, "Set3"), brewer.pal(nrow(chromosome_data)-12, "Dark2"))
    names(nucl_col) <- as.character(chromosome_data$Chromosome)
} else if(nrow(chromosome_data)>20 && nrow(chromosome_data)<=23) {
    nucl_col <- c(brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"), brewer.pal(nrow(chromosome_data)-20, "Set1"))
    names(nucl_col) <- as.character(chromosome_data$Chromosome)
}

length_sector <- as.numeric(length(sector))

svg(out_svg, width = 9, height = 7, bg = "transparent")

circos.par(cell.padding = c(0, 0, 0, 0), start.degree = 270, gap.degree = c(rep(1,length_sector-3),10, 1, 10))
circos.initialize(factors = factor(sector, levels = sector), xlim = sector_xlim,
sector.width = c(sector_xlim[1:(length_sector-2),2]/sum(sector_xlim[1:(length_sector-2),2]),
                    0.5*sector_xlim[(length_sector-1),2]/sum(sector_xlim[(length_sector-1),2]),
                    0.5*sector_xlim[length_sector,2]/sum(sector_xlim[length_sector,2])))

circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(mean(xlim), mean(ylim),
                sector.index,
                cex = 0.8, font = 2,
                facing = "bending.inside",
                niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA)

circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    if(sector.index=="MT"){
        circos.axis(h = "top", labels.cex = 0.5)
        circos.par("track.margin")
    } else {
        circos.axis(h = "top", labels = FALSE, major.tick = FALSE, minor.ticks = FALSE)
        circos.par("track.margin" = c(0, 0.01))
    }
    
}, track.height = 0.02, bg.col = c(nucl_col, mt_col), track.margin = c(0, 0.02))

circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
}, track.height = 0.02, track.margin = c(0, 0.01))

for(i in 1:nrow(final_data)) {
    circos.link(final_data[i,4], c(final_data[i,5], final_data[i,6]),
            final_data[i,1], c(final_data[i,2], final_data[i,3]),
            col = paste0(nucl_col[final_data[i,1]], "30"), border = "red", lwd = 0.2)
    circos.rect(final_data[i,3], 0, final_data[i,2], 1, sector.index = final_data[i,1], col = mt_col[final_data[i,4]])
    circos.rect(final_data[i,6], 0, final_data[i,5], 1, sector.index = final_data[i,4], col = nucl_col[final_data[i,1]])
    
}

dev.off()

circos.clear()

png(out_png, res = 600, width = 9, height = 6, units = "in", bg="transparent")

circos.par(cell.padding = c(0, 0, 0, 0), start.degree = 270, gap.degree = c(rep(1,length_sector-3),10, 1, 10))
circos.initialize(factors = factor(sector, levels = sector), xlim = sector_xlim,
sector.width = c(sector_xlim[1:(length_sector-2),2]/sum(sector_xlim[1:(length_sector-2),2]),
                    0.5*sector_xlim[(length_sector-1),2]/sum(sector_xlim[(length_sector-1),2]),
                    0.5*sector_xlim[length_sector,2]/sum(sector_xlim[length_sector,2])))

circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(mean(xlim), mean(ylim),
                sector.index,
                cex = 0.8, font = 2,
                facing = "bending.inside",
                niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA)

circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    if(sector.index=="MT"){
        circos.axis(h = "top", labels.cex = 0.5)
        circos.par("track.margin")
    } else {
        circos.axis(h = "top", labels = FALSE, major.tick = FALSE, minor.ticks = FALSE)
        circos.par("track.margin" = c(0, 0.01))
    }
    
}, track.height = 0.02, bg.col = c(nucl_col, mt_col), track.margin = c(0, 0.02))

circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
}, track.height = 0.02, track.margin = c(0, 0.01))

for(i in 1:nrow(final_data)) {
    circos.link(final_data[i,4], c(final_data[i,5], final_data[i,6]),
            final_data[i,1], c(final_data[i,2], final_data[i,3]),
            col = paste0(nucl_col[final_data[i,1]], "30"), border = "red", lwd = 0.2)
    circos.rect(final_data[i,3], 0, final_data[i,2], 1, sector.index = final_data[i,1], col = mt_col[final_data[i,4]])
    circos.rect(final_data[i,6], 0, final_data[i,5], 1, sector.index = final_data[i,4], col = nucl_col[final_data[i,1]])
    
}

dev.off()

circos.clear()
