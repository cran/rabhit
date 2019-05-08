# HaplotypeR graphic functions -----------------------------------------------------

#' @include rabhit.R
#' @include internal_functions.R
NULL

#' Graphical output of an inferred haplotype
#'
#' The \code{plotHaplotype} functions visualizes an inferred haplotype.
#'
#'
#' @param    hap_table            haplotype summary table. See details.
#' @param    html_output          if TRUE, a html5 interactive graph is outputed. Defualt is FALSE.
#' @param    gene_sort            if by 'name' the genes in the output are ordered lexicographically,
#' if by 'position' only functional genes are used and are ordered by their chromosomal location. Default is 'position'.
#' @param    text_size            the size of graph labels. Default is 14 (pts).
#' @param    removeIGH            if TRUE, 'IGH'\'IGK'\'IGL' prefix is removed from gene names.
#' @param    plotYaxis            if TRUE, Y axis labels (gene names) are plotted on the middle and right plots. Default is TRUE.
#' @param    chain                the Ig chain: IGH,IGK,IGL. Default is IGH.
#' @param    dir                  The output folder for saving the haplotype map for multiple individuals.
#'
#' @return
#'
#' A haplotype map visualization. If more than one subject is visualized, a pdf is created. If html_output is TRUE, a folder named html_output is created with individual graphs.
#'
#' @details
#'
#' A \code{data.frame} in a haplotype format created by \code{createFullHaplotype} function.
#'
#' @examples
#'
#' # Selecting a single individual from the haplotype samples data
#' haplo_db = samplesHaplotype[samplesHaplotype$SUBJECT=='I5', ]
#'
#' # plot haplotype
#' plotHaplotype(haplo_db)
#'
#' @export
plotHaplotype <- function(hap_table, html_output = FALSE, gene_sort = c("name", "position"), text_size = 14, removeIGH = TRUE, plotYaxis = TRUE, chain = c("IGH",
    "IGK", "IGL"), dir) {
    if (missing(chain)) {
        chain = "IGH"
    }
    chain <- match.arg(chain)

    if (missing(gene_sort)) {
        gene_sort = "position"
    }
    gene_sort <- match.arg(gene_sort)

    hapBy_cols = names(hap_table)[grep(chain, names(hap_table))]

    hapBy_alleles = gsub("_", "*", hapBy_cols)

    if (!("SUBJECT" %in% names(hap_table))) {
        hap_table$SUBJECT <- rep("S1", nrow(hap_table))
    }

    plot_list <- c()
    for (sample_name in unique(hap_table$SUBJECT)) {


        GENE.loc.tmp <- GENE.loc[[chain]]

        haplo.db <- parseHapTab(hap_table[hap_table$SUBJECT == sample_name, ], chain = chain)
        geno.df <- sortDFByGene(haplo.db$geno.df, chain = chain, method = gene_sort, removeIGH = removeIGH)
        kval.df <- sortDFByGene(haplo.db$kval.df, chain = chain, method = gene_sort, removeIGH = removeIGH)
        count.df <- sortDFByGene(haplo.db$count.df, chain = chain, method = gene_sort, removeIGH = removeIGH)

        ########################################################################################################

        ### Prepare All panels

        if (length(grep("[0-9][0-9]_[0-9][0-9]", geno.df$ALLELES)) != 0) {
            geno.df <- as.data.frame(geno.df %>% group_by(.data$hapBy, .data$GENE) %>% mutate(n = n()))
            geno.df$freq <- ifelse(geno.df$n == 2, 0.5, ifelse(geno.df$n != 1, 0.25, 1))
            non_reliable_alleles_text <- nonReliableAllelesText(geno.df[grep("[0-9][0-9]_[0-9][0-9]", geno.df$ALLELES), ])
            geno.df$ALLELES[grep("[0-9][0-9]_[0-9][0-9]", geno.df$ALLELES)] <- "NRA"

            allele_palette <- alleleHapPalette(geno.df$ALLELES)
            AlleleCol <- allele_palette$AlleleCol
            transper <- allele_palette$transper
        } else {
            allele_palette <- alleleHapPalette(geno.df$ALLELES, NRA = F)
            AlleleCol <- allele_palette$AlleleCol
            transper <- allele_palette$transper
            non_reliable_alleles_text <- c()
        }

        ## Middle panel
        geno.df$ALLELES <- factor(geno.df$ALLELES, levels = AlleleCol)
        p = ggplot(geno.df, aes_string(x = "GENE", fill = "ALLELES")) + geom_bar(position = "fill", width = 0.9, na.rm = T) + coord_flip() +
            xlab("") + ylab("") + facet_grid(paste0(".~", "hapBy"), switch = "x") + theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = text_size), strip.background = element_blank(),
            strip.text = element_text(face = "bold"), axis.text = element_text(colour = "black"), panel.spacing = unit(0, "cm"), strip.switch.pad.grid = unit(0,
                "cm"), plot.margin = unit(c(0.25, 0, 0.2, 0), "cm"), legend.key = element_rect("#DCDCDC")) + scale_fill_manual(values = alpha(names(AlleleCol), transper), name = "Alleles", drop = F)

        if (!plotYaxis) {
            p = p + theme(axis.text.y = element_blank())
        }

        ## Right panel plot K values

        pk <- ggplot(kval.df, aes_string(x = "GENE", fill = "K_GROUPED")) + theme_bw() + theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = text_size), strip.background = element_blank(),
            strip.text = element_text(face = "bold"), panel.spacing = unit(0, "cm"), strip.switch.pad.grid = unit(0, "cm"), plot.margin = unit(c(0.25, 0,
                0.2, 0), "cm"), legend.key = element_rect("#DCDCDC")) + geom_bar(position = "fill", width = 0.7, na.rm = T) + coord_flip() + xlab("") + ylab("") + facet_grid(paste0(".~", "hapBy"),
            switch = "x")

        count.df$ALLELES <- as.character(count.df$ALLELES)
        count.df$ALLELES[grep("[0-9][0-9]_[0-9][0-9]", count.df$ALLELES)] <- "NRA"
        count.df$border <- factor(ifelse(count.df$ALLELES == "NRA", "black", "white"), levels = c("black", "white"))
        count.df$ALLELES <- factor(count.df$ALLELES, levels = AlleleCol)
        ## Left panel
        p2 <- ggplot(count.df, aes_string(x = "GENE", y = "COUNT2", fill = "ALLELES")) + geom_bar(stat = "identity", position = "Dodge", width = 0.9,
            na.rm = T, aes_string(colour = "border")) + coord_flip() + background_grid(minor = "none") + scale_fill_manual(values = alpha(names(AlleleCol),
            transper), name = "ALLELES", drop = F) + scale_color_manual(values = alpha(c("black", "white"), c(0.5, 0)), drop = F) + theme(legend.position = "none",
            strip.text = element_text(face = "bold"), axis.text = element_text(colour = "black"), text = element_text(size = text_size), plot.margin = unit(c(0.25,
                0, -0.05, 0), "cm"), panel.background = element_blank(), legend.key = element_rect("#DCDCDC")) + scale_y_continuous(breaks = seq(-3, 3, by = 1), labels = c(3:0, 1:3)) + ylab(expression("log"[10] *
            "(Count+1)")) + xlab("Gene") + geom_hline(yintercept = c(0), linetype = "dotted")

        if (is.data.frame(non_reliable_alleles_text)) {
            p <- p + geom_text(data = non_reliable_alleles_text, aes_string(label = "text", x = "GENE", y = "pos"), angle = 0, size = non_reliable_alleles_text$size)
            non_reliable_alleles_text$COUNT2 <- count.df$COUNT2[count.df$ALLELES == "NRA"]
            non_reliable_alleles_text <- non_reliable_alleles_text[non_reliable_alleles_text$COUNT2 != 0, ]
            non_reliable_alleles_text$hjust <- ifelse(non_reliable_alleles_text$COUNT2 >= 0, 0, 1)
            p2 <- p2 + geom_text(data = non_reliable_alleles_text, aes_string(label = "text", x = "GENE", hjust = "hjust"), angle = 0, size = 2.25)
        }
        ########################################################################################################

        ### Plot All panels

        if (html_output) {

            ## Prepare panels for html plot

            p = p + theme(axis.title.x = element_blank())
            p.l <- ggplotly(p, height = 800, width = 400) %>% plotly::layout(showlegend = FALSE)

            pk = pk + theme(axis.title = element_blank())
            pk <- pk + scale_fill_manual(name = "log<sub>10</sub>(lK)", values = c('#FFFFFF',RColorBrewer::brewer.pal(9,'Blues')), drop = FALSE)
            pk.l <- ggplotly(pk, height = 800, width = 400) %>% plotly::layout(showlegend = TRUE)
            pk.l$x$layout$annotations[[1]]$text = p.l$x$layout$annotations[[1]]$text
            pk.l$x$layout$annotations[[2]]$text = p.l$x$layout$annotations[[2]]$text

            p2 <- p2 + ylab("log<sub>10</sub>(Count+1)")
            p2.l <- ggplotly(p2, height = 1000, width = 700) %>% plotly::layout(margin = list(b = 50), yaxis = list(title = paste0(c(rep("&nbsp;", 3), "Gene",
                rep("&nbsp;", 3), rep("\n&nbsp;", 1)), collapse = "")), showlegend = TRUE)

            p2.l$x$layout$xaxis$ticktext = c(lapply(p2.l$x$layout$xaxis$ticktext[1:match("0", p2.l$x$layout$xaxis$ticktext) - 1], function(x) paste0("-",
                x)), p2.l$x$layout$xaxis$ticktext[match("0", p2.l$x$layout$xaxis$ticktext):length(p2.l$x$layout$xaxis$ticktext)])

            p2.l$x$data[[length(p2.l$x$data)]]$y[2] = p2.l$x$layout$yaxis$range[2]

            add_border <- function(p_l){
              p_l$x$attrs <- lapply(p_l$x$attrs,
                                      function(x){
                                        x$mode <- "markers"
                                        x
                                      })

              for (i in 1:length(p_l$x$data)) {
                p_l$x$data[[i]]$marker$line$width <- suppressMessages(0.2)
                p_l$x$data[[i]]$marker$line$color <- ifelse(p_l$x$data[[i]]$marker$line$color=='transparent', 'black', p_l$x$data[[i]]$marker$line$color)
              }
              return(p_l)
            }

            p.l <- add_border(p.l)
            pk.l <- add_border(pk.l)
            mgsub <- function(pattern, replacement, x, ...) {
                if (length(pattern) != length(replacement)) {
                  stop("pattern and replacement do not have the same length.")
                }
                result <- x
                for (i in 1:length(pattern)) {
                  result <- gsub(pattern[i], replacement[i], result, ...)
                }
                result
            }


            text_for_hovertext <- function(labels, count.df) {
                for (i in 1:length(labels)) {
                  label <- labels[i]
                  gene <- strsplit(strsplit(label, "<")[[1]][1], " ")[[1]][2]
                  allele <- strsplit(label, "Allele: ")[[1]][2]
                  if (!is.na(NA)) {
                    count <- strsplit(strsplit(label, "<br />Count: ")[[1]][2], "<")[[1]][1]
                    if (count == "NA")
                      next
                    count <- as.numeric(count)
                    if (count%%1 != 0)
                      count <- count.df %>% filter(.data$GENE == gene & .data$ALLELES == allele & round(.data$COUNT3, nchar(as.character(count)) - 2) == count) %>% select(.data$COUNT) else count <- count.df %>% filter(.data$GENE == gene & .data$ALLELES == allele & .data$COUNT3 == count) %>% select(.data$COUNT)
                    labels[i] <- paste0("Gene: ", gene, "<br />Allele: ", allele, "<br />Count: ", count[1, ])
                  }
                }
                return(labels)
            }

            text_for_hovertext_non_reliable <- function(labels, text) {
              for (i in 1:length(labels)) {
                label <- labels[i]
                gene <- strsplit(label, "GENE: ")[[1]][2]
                allele <- gsub("text: ","",strsplit(grep(gene,text,value=T), "<")[[1]][1])
                labels[i] <- paste0("Gene: ", gene, "<br />Allele: ", allele)
              }
              return(labels)
            }

            count.df$COUNT3 <- abs(count.df$COUNT2)

            for (i in 1:length(p2.l$x$data)) {
                p2.l$x$data[[i]]$text <- mgsub(c("~GENE", "~COUNT2", "~factor[(]ALLELES[,] levels [=] AlleleCol[)]", "~yintercept: 0"), c("Gene", "Count",
                  "Allele", ""), p2.l$x$data[[i]]$text)
                p2.l$x$data[[i]]$text <- text_for_hovertext(p2.l$x$data[[i]]$text, count.df)
            }

            if(is.data.frame(non_reliable_alleles_text)){
              ind <- grep('NRA',p.l$x$data)
              print(length(ind))
              if(length(ind) == 4){
                p.l$x$data[[ind[1]]]$text <- text_for_hovertext_non_reliable(p.l$x$data[[ind[1]]]$text, p.l$x$data[[ind[3]]]$hovertext)
                p.l$x$data[[ind[2]]]$text <- text_for_hovertext_non_reliable(p.l$x$data[[ind[2]]]$text, p.l$x$data[[ind[4]]]$hovertext)
                p.l$x$data[[ind[3]]] <- NULL
                p.l$x$data[[ind[4]-1]] <- NULL
              }else{
                p.l$x$data[[ind[1]]]$text <- text_for_hovertext_non_reliable(p.l$x$data[[ind[1]]]$text, p.l$x$data[[ind[2]]]$hovertext)
                p.l$x$data[[ind[2]]] <- NULL
              }
            }

            for (i in 1:length(p.l$x$data)) {
                p.l$x$data[[i]]$text <- mgsub(c("~GENE", "~factor[(]ALLELES[,] levels [=] AlleleCol[)]"), c("Gene", "Allele"), p.l$x$data[[i]]$text)
            }

            for (i in 1:length(pk.l$x$data)) {
                pk.l$x$data[[i]]$text <- mgsub(c("~GENE", "~K[_]GROUPED"), c("Gene", "K"), pk.l$x$data[[i]]$text)
            }


            p.l.c <- suppressWarnings(subplot(p2.l, p.l, pk.l, widths = c(0.4, 0.2, 0.2), shareY = T, titleX = TRUE, margin = 0.01, which_layout = 1))


            for (i in 1:length(p2.l$x$data)) {
                p.l.c$x$data[[i]]$showlegend <- FALSE
            }

            for(i in 1:(length(p2.l$x$data)+length(p.l$x$data))){
              p.l.c$x$data[[i]]$name <- gsub("[()|,]|white|black","",p.l.c$x$data[[i]]$name)
            }

            p.l.c$x$layout$annotations[[6]]$text = "log<sub>10</sub>(lK)"
            p.l.c$x$layout$annotations[[6]]$xanchor = "center"
            p.l.c$x$layout$annotations[[6]]$y = 1 - 0.035 * (length(grep('[[]',grep('TRUE',p.l.c$x$data,value = T),invert = T))) #(length(AlleleCol) + 2.4)  #0.52
            p.l.c$x$layout$annotations[[6]]$x = 0.98
            p.l.c$x$layout$annotations[[6]]$font$size = 16


            p.l.c$x$layout$annotations[[3]] <- p.l.c$x$layout$annotations[[6]]
            p.l.c$x$layout$annotations[[3]]$text = "Alleles"
            p.l.c$x$layout$annotations[[3]]$y = 0.99
            p.l.c$x$layout$annotations[[3]]$x = 0.98
            p.l.c$x$layout$annotations[[3]]$legendTitle = FALSE
            p.l.c$x$layout$annotations[[3]]$font$size = 16

            plot_list[[sample_name]] <- p.l.c

        } else {
            p.legend <- get_legend(p + theme(legend.key = element_rect("#DCDCDC")))
            p = p + theme(legend.position = "none",axis.title.x = element_blank())

            pk = pk + scale_fill_manual(name = expression("log"[10] * "(lK)"), values = c('#FFFFFF',RColorBrewer::brewer.pal(9,'Blues')), drop = FALSE)
            pk.legend <- get_legend(pk + theme(legend.key = element_rect("gray")))
            pk = pk + theme(legend.position = "none",axis.title = element_blank())

            p.legends <- plot_grid(pk.legend, p.legend, ncol = 1, align = "hv")

            p1 <- plot_grid(p2, p, pk, nrow = 1, rel_widths = c(0.35, 0.15, 0.05), align = "hv", axis = "b")

            p <- plot_grid(p1, p.legends, ncol = 2, rel_widths = c(1, 0.1))

            plot_list[[sample_name]] <- p
        }

    }
    if (length(plot_list) != 1) {


        if(!missing(dir)){
          dir <- file.path(dir, "haplotype_output")
          dir.create(dir)

        }else{
          dir <- tempdir()
        }

        if (html_output) {

            for (sample_name in names(plot_list)) {
                htmlwidgets::saveWidget(plot_list[[sample_name]], paste0(dir, '/', sample_name, ".html"), selfcontained = T)
            }
        } else {
            pdf(paste0(dir, "/haplotype_output.pdf"), height = 20, width = 15)
            for (sample_name in names(plot_list)) {

                title <- ggdraw() + draw_label(sample_name, fontface = "bold")
                plot(plot_grid(title, plot_list[[sample_name]], ncol = 1, rel_heights = c(0.05, 1)))
            }

            dev.off()
        }

    } else if (html_output)
        return(plot_list[[1]]) else plot(plot_list[[1]])
}

########################################################################################################
#' Graphical output of alleles division by chromosome
#'
#' The \code{hapHeatmap} function generates a graphical output of the alleles per gene in multiple samples.
#'
#'
#' @param    hap_table            haplotype summary table. See details.
#' @param    chain                the IG chain: IGH,IGK,IGL. Default is IGH.
#' @param    gene_sort            if by 'name' the genes in the output are ordered lexicographically,
#' if by 'position' only functional genes are used and are ordered by their chromosomal location. Default is 'position'.
#' @param    removeIGH            if TRUE, 'IGH'\'IGK'\'IGL' prefix is removed from gene names.
#' @param    lk_cutoff            the lK cutoff value to be considerd low for texture layer. Defualt is lK<1.
#' @param    mark_low_lk          if TRUE, a texture is add for low lK values. Defualt is TRUE.
#'
#' @return
#'
#' A heat-map visualization of the haplotype inference for multiple samples.
#'
#' @details
#'
#' A \code{data.frame} created by \code{createFullHaplotype}.
#'
#' @examples
#' # Plotting haplotpe heatmap
#' hapHeatmap(samplesHaplotype)
#'
#' @export
hapHeatmap <- function(hap_table, chain = c("IGH", "IGK", "IGL"), gene_sort = "position", removeIGH = TRUE, lk_cutoff = 1, mark_low_lk = TRUE) {


    if (missing(chain)) {
        chain = "IGH"
    }
    chain <- match.arg(chain)

    hapBy_alleles <- gsub("_", "*", names(hap_table)[c(3,4)])
    hapBy_cols <- gsub("IG[H|K|L]", "", hapBy_alleles)
    samples <- unique(hap_table$SUBJECT)

    haplo_db <- c()
    for (sample_name in samples) {
      haplo.db <- parseHapTab(hap_table[hap_table$SUBJECT == sample_name, ], chain = chain)
      geno.df <- sortDFByGene(haplo.db$geno.df, chain = chain, method = gene_sort, removeIGH = removeIGH)
      kval.df <- sortDFByGene(haplo.db$kval.df, chain = chain, method = gene_sort, removeIGH = removeIGH)
      geno.df$K <- apply(geno.df[, c("GENE", "hapBy")], 1, function(x) {
        asNum(kval.df$K[kval.df$GENE == x[[1]] & kval.df$hapBy == x[[2]]], na.strings = "NA")})
      haplo_db <- rbind(haplo_db, geno.df)
    }


    haplo_db_texture <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), c('SUBJECT','GENE','ALLELES','hapBy','K','points','yend','x','xend'))
    loc <- 1:length(levels(droplevels(haplo_db$GENE)))
    names(loc) <- rev(levels(droplevels(haplo_db$GENE)))
    for (i in 1:nrow(haplo_db)) {
        if (haplo_db$K[i] < lk_cutoff && !haplo_db$ALLELES[i] %in% c("Unk", "Del", "NR")) {
            tmp_point <- haplo_db[i, ] %>% slice(rep(1, each = ifelse(length(samples) < 4, 15, 8))) %>% mutate(points = seq(0, 0.9, length.out = ifelse(length(samples) <
                4, 15, 8)), yend = seq(0, 0.9, length.out = ifelse(length(samples) < 4, 15, 8)) + 0.1, x = loc[as.character(.data$GENE)] - 0.49, xend = loc[as.character(.data$GENE)] + 0.49)
            haplo_db_texture <- bind_rows(haplo_db_texture, tmp_point)
        }
    }


    heatmap.df <- as.data.frame(haplo_db %>% group_by(.data$SUBJECT, .data$hapBy, .data$GENE) %>% mutate(n = n()))
    heatmap.df$freq <- ifelse(heatmap.df$n == 2, 0.5, ifelse(heatmap.df$n != 1, 0.25, 1))
    heatmap.df$GENE <- factor(heatmap.df$GENE, levels = gsub("IG[H|K|L]", "", GENE.loc[[chain]]))
    heatmap.df$title <- ifelse(heatmap.df$hapBy == hapBy_cols[1], hapBy_cols[1], hapBy_cols[2])


    if (length(grep("[0-9][0-9]_[0-9][0-9]", heatmap.df$ALLELES)) != 0) {
        non_reliable_alleles_text <- nonReliableAllelesText(heatmap.df[grep("[0-9][0-9]_[0-9][0-9]", heatmap.df$ALLELES), ])
        heatmap.df$ALLELES[grep("[0-9][0-9]_[0-9][0-9]", heatmap.df$ALLELES)] <- "NRA"
        allele_palette <- alleleHapPalette(heatmap.df$ALLELES)
        non_reliable_alleles_text$ALLELES <- factor(non_reliable_alleles_text$ALLELES, levels = allele_palette$AlleleCol)
    } else {
        non_reliable_alleles_text <- c()
        allele_palette <- alleleHapPalette(heatmap.df$ALLELES, NRA = F)

    }
    sub_lev <- F
    if (length(levels(hap_table$SUBJECT)) != 0) {
        sub_lev <- T
        heatmap.df$SUBJECT <- factor(heatmap.df$SUBJECT, levels = levels(hap_table$SUBJECT))
        non_reliable_alleles_text$SUBJECT <- factor(non_reliable_alleles_text$SUBJECT, levels = levels(hap_table$SUBJECT))
    }

    heatmap.df$ALLELES <- factor(heatmap.df$ALLELES, levels = allele_palette$AlleleCol)
    col <- ifelse(c(1:(length(unique(heatmap.df[heatmap.df$hapBy == hapBy_cols[1], "SUBJECT"])) * 4))%%3 == 0, "black", "transparent")
    p <- ggplot() + geom_col(data = heatmap.df[heatmap.df$hapBy == hapBy_cols[1], ], mapping = aes_string(x = "GENE", y = "freq", fill = "ALLELES"), position = "fill", width = 0.95, na.rm = T) +
        scale_fill_manual(values = alpha(names(allele_palette$AlleleCol), allele_palette$transper), name = "Alleles", drop = FALSE) + facet_grid(SUBJECT ~
        title, as.table = FALSE, switch = "y") + scale_y_continuous(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) + theme(axis.text.x = element_text(angle = 90,
        vjust = 0.5, hjust = 1, size = 10, colour = "black"), strip.text.x = element_text(size = 14), strip.text.y = element_text(angle = 180, size = 10),
        strip.placement = "outside", strip.background =element_rect(fill="seashell2"), axis.ticks.y = element_line(colour = col), axis.text.y = element_blank(), strip.background.y = element_blank(), legend.position = "bottom",
        legend.justification = "center", panel.spacing.y = unit(0.9, "pt")) + labs(y = "", x = "") + guides(fill = guide_legend(nrow = round(length(allele_palette$AlleleCol)/10),
        order = 1, override.aes = list(color = "#DCDCDC")))

    col <- ifelse(c(1:(length(unique(heatmap.df[heatmap.df$hapBy == hapBy_cols[2], "SUBJECT"])) * 4))%%3 == 0, "black", "transparent")
    p1 <- ggplot() + geom_col(data = heatmap.df[heatmap.df$hapBy == hapBy_cols[2], ], mapping = aes_string(x = "GENE", y = "freq", fill = "ALLELES"), position = "fill",
        width = 0.95, na.rm = T) + scale_fill_manual(values = alpha(names(allele_palette$AlleleCol), allele_palette$transper), name = "Alleles", drop = FALSE) +
        facet_grid(SUBJECT ~ title, as.table = FALSE, switch = "y") + scale_y_continuous(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) + theme(axis.text.x = element_text(angle = 90,
        vjust = 0.5, hjust = 1, size = 10, colour = "black"), strip.text.x = element_text(size = 14), strip.text.y = element_text(angle = 180, size = 10),
        strip.placement = "outside", strip.background =element_rect(fill="seashell2"), axis.ticks.y = element_line(colour = col), axis.text.y = element_blank(), strip.background.y = element_blank(), legend.position = "bottom",
        legend.justification = "center", panel.spacing.y = unit(0.9, "pt")) + labs(y = "") + guides(fill = guide_legend(nrow = round(length(allele_palette$AlleleCol)/10),
        order = 1, override.aes = list(color = "#DCDCDC")))


    if (mark_low_lk & nrow(haplo_db_texture) != 0) {

        haplo_db_texture <- haplo_db_texture[!duplicated(haplo_db_texture[, c("GENE", "ALLELES", "hapBy", "K", "points", "SUBJECT")]), ]
        haplo_db_texture$col <- "<1"
        if(sub_lev) haplo_db_texture$SUBJECT <- factor(haplo_db_texture$SUBJECT, levels = levels(hap_table$SUBJECT))
        haplo_db_texture$ALLELES <- factor(haplo_db_texture$ALLELES, levels = allele_palette$AlleleCol)
        haplo_db_texture$ALLELES[grep("[0-9][0-9]_[0-9][0-9]", haplo_db_texture$ALLELES)] <- "NRA"

        if (gsub(chain,"",hapBy_alleles[1]) %in% haplo_db_texture$hapBy) {
            # Get Allele legend
            gt1 = ggplotGrob(p)

            p <- p + geom_segment(data = haplo_db_texture[haplo_db_texture$hapBy == hapBy_cols[1], ], mapping = aes_string(x = "x", xend = "xend", y = "points", yend = "yend", color = "col"), colour = "white")

            if (gsub(chain,"",hapBy_alleles[2]) %in% haplo_db_texture$hapBy) {
                p1 <- p1 + geom_segment(data = haplo_db_texture[haplo_db_texture$hapBy == hapBy_cols[2], ], mapping = aes_string(x = "x", xend = "xend", y = "points", yend = "yend", color = "col") , colour = "white")
            }

            gt2 = getDigLegend(unique(haplo_db_texture$col))

        } else {
            # Get Allele legend
            gt1 = ggplotGrob(p1)

            p1 <- p1 + geom_segment(data = haplo_db_texture[haplo_db_texture$hapBy == hapBy_cols[2], ], mapping = aes_string(x = "x", xend = "xend", y = "points", yend = "yend", color = "col"), colour = "white")
            gt2 = getDigLegend(unique(haplo_db_texture$col))
        }


        # Get lK legend
        leg1 = gtable::gtable_filter(gt1, "guide-box")
        leg2 = gtable::gtable_filter(gt2, "guide-box")
        # Combine the legends
        leg <- cbind(leg1[["grobs"]][[1]], leg2[["grobs"]][[1]], size = "first")
        # Insert legend into g1 (or g2)
        gt1$grobs[gt1$layout$name == "guide-box"][[1]] <- leg
        gt1$grobs[gt1$layout$name == "guide-box"][[1]]$layout[3, c("t", "b")] <- gt1$grobs[gt1$layout$name == "guide-box"][[1]]$layout[1, c("t", "b")]
        gt1$grobs[gt1$layout$name == "guide-box"][[1]]$layout[3, c("l", "r")] <- gt1$grobs[gt1$layout$name == "guide-box"][[1]]$layout[1, c("l", "r")] +
            2
        legend <- get_legend(gt1)
    } else legend <- get_legend(p)

    if (is.data.frame(non_reliable_alleles_text)) {
      p <- p + geom_text(data = non_reliable_alleles_text[non_reliable_alleles_text$hapBy == hapBy_cols[1], ], aes_string(label = "text", x = "GENE", y = "pos"), angle = 90,
                         size = non_reliable_alleles_text$size[non_reliable_alleles_text$hapBy == hapBy_cols[1]])
      p1 <- p1 + geom_text(data = non_reliable_alleles_text[non_reliable_alleles_text$hapBy == hapBy_cols[2], ], aes_string(label = "text", x = "GENE", y = "pos"),
                           angle = 90, size = non_reliable_alleles_text$size[non_reliable_alleles_text$hapBy == hapBy_cols[2]])
    }

    plot(plot_grid(plot_grid(p + theme(legend.position = "none", plot.background = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill = "transparent", colour = NA), axis.line = element_line(color="black")),
                             p1 + theme(legend.position = "none", plot.background = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill = "transparent", colour = NA),axis.line = element_line(color="black")), nrow = 2, align = "hv"), legend, rel_heights = c(0.9,
        0.1), ncol = 1))
}

########################################################################################################
#' Hierarchical clustering of haplotypes graphical output
#'
#' The \code{hapDendo} function generates a graphical output of an hierarchical clustering based on the Jaccard distance between multiple samples' haplotypes.
#'
#'
#' @param    hap_table            haplotype summary table. See details.
#' @param    chain                the IG chain: IGH,IGK,IGL. Default is IGH.
#' @param    gene_sort            if by 'name' the genes in the output are ordered lexicographically,
#' if by 'position' only functional genes are used and are ordered by their chromosomal location. Default is 'position'.
#' @param    removeIGH            if TRUE, 'IGH'\'IGK'\'IGL' prefix is removed from gene names. Defualt is TRUE.
#' @param    mark_low_lk          if TRUE, a texture is add for low lK values. Defualt is TRUE.
#' @param    lk_cutoff            the lK cutoff value to be considerd low for texture layer. Defualt is lK<1.
#'
#' @return
#'
#' A multitple samples visualization of the distances between haplotypes.
#'
#' @details
#'
#' A \code{data.frame} created by \code{createFullHaplotype}.
#'
#' @examples
#' # Plotting haplotype hierarchical clustering based on the Jaccard distance
#' hapDendo(samplesHaplotype)
#'
#' @export
hapDendo <- function(hap_table, chain = c("IGH", "IGK", "IGL"), gene_sort = c("name", "position"), removeIGH = TRUE, mark_low_lk = TRUE, lk_cutoff = 1) {

  if (missing(chain)) {
    chain = "IGH"
  }
  chain <- match.arg(chain)

  if (missing(gene_sort)) {
    gene_sort = "position"
  }
  gene_sort <- match.arg(gene_sort)
  hapBy_cols = names(hap_table)[c(3,4)]
  samples <- unique(hap_table$SUBJECT)
  if (length(samples) < 2)
    stop("hapDendo function requires at least two samples")
  # creating the distance matrix for clustering
  mat <- matrix(NA, nrow = length(samples), ncol = length(samples))
  for (i in 2:length(samples)) {
    for (j in 1:(i - 1)) {

      hap_merge <- merge(hap_table[hap_table$SUBJECT == samples[i], c("GENE", hapBy_cols[1], hapBy_cols[2])], hap_table[hap_table$SUBJECT == samples[j],
                                                                                                                        c("GENE", hapBy_cols[1], hapBy_cols[2])], by = "GENE")

      mat[i, j] <- calcJacc(vec1A = hap_merge[, 2], vec1B = hap_merge[, 3], vec2A = hap_merge[, 4], vec2B = hap_merge[, 5], method = "geneByGene")

    }
  }
  colnames(mat) <- samples
  rownames(mat) <- samples

  # finding the hierarchical clustering
  fit <- hclust(as.dist(mat), method = "ward.D")
  samples <- samples[fit$order]

  # Preparing the hclust data for plotting
  dend <- as.dendrogram(hclust(as.dist(mat), method = "ward.D"))
  dend_data <- ggdendro::dendro_data(dend)
  segment_data <- with(ggdendro::segment(dend_data), data.frame(x = y, y = x, xend = yend, yend = xend))

  # Using the dendrogram label data to position the samples labels
  samples_pos_table <- with(dend_data$labels, data.frame(y_center = x, gene = as.character(label), height = 1))
  samples_axis_limits <- with(samples_pos_table, c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))) + 0.1 * c(-1, 1)

  plt_dendr <- ggplot(segment_data) + geom_segment(aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) + scale_x_continuous(expand = c(0, 0.01)) + scale_y_continuous(breaks = samples_pos_table$y_center,
                                                                                                                                                                              labels = samples_pos_table$gene, limits = samples_axis_limits, expand = c(0, 0)) + labs(x = "Jaccard distance", y = "", colour = "", size = "") +
    theme_bw() + theme(panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),
                       axis.text = element_text(size = 14, colour = "black"), axis.title.x = element_text(size = 14, colour = "black")) + ylab('')

  # Creating the haploype db for ploting
  haplo_db_clust <- c()
  for (sample_name in samples) {
    haplo.db <- parseHapTab(hap_table[hap_table$SUBJECT == sample_name, ], chain = chain)
    geno.df <- sortDFByGene(haplo.db$geno.df, chain = chain, method = gene_sort, removeIGH = removeIGH)
    kval.df <- sortDFByGene(haplo.db$kval.df, chain = chain, method = gene_sort, removeIGH = removeIGH)
    geno.df$K <- apply(geno.df[, c("GENE", "hapBy")], 1, function(x) {
      asNum(kval.df$K[kval.df$GENE == x[[1]] & kval.df$hapBy == x[[2]]], na.strings = "NA")})
    haplo_db_clust <- rbind(haplo_db_clust, geno.df)
  }
  allele_palette <- alleleHapPalette(haplo_db_clust$ALLELES)

  # Formating the data to fit heatmap plot
  haplo_db_clust <- haplo_db_clust %>% group_by(.data$SUBJECT, .data$hapBy, .data$GENE) %>% mutate(n = n())
  haplo_db_clust$freq <- ifelse(haplo_db_clust$n == 2, 0.5, ifelse(haplo_db_clust$n != 1, 0.25, 1))
  haplo_db_clust$GENE <- factor(haplo_db_clust$GENE, levels = gsub("IG[H|K|L]", "", GENE.loc[[chain]]))
  haplo_db_clust$grouper_x <- "Gene"
  haplo_db_clust$grouper_y <- sapply(1:nrow(haplo_db_clust), function(i) paste0(haplo_db_clust$SUBJECT[i], " ", haplo_db_clust$hapBy[i]))

  haplo_db_clust_texture <- setNames(data.frame(matrix(ncol = 13, nrow = 0)), c('SUBJECT','GENE','ALLELES','hapBy','K','n','freq','grouper_x','grouper_y','points','yend','x','xend'))
  loc <- 1:length(levels(droplevels(haplo_db_clust$GENE)))
  names(loc) <- levels(droplevels(haplo_db_clust$GENE))
  for (i in 1:nrow(haplo_db_clust)) {
    if (haplo_db_clust$K[i] < lk_cutoff && !haplo_db_clust$ALLELES[i] %in% c("Unk", "Del", "NR")) {
      tmp_point <- haplo_db_clust[i, ] %>% slice(rep(1, each = ifelse(length(samples) < 4, 15, 8))) %>% mutate(points = seq(0, 0.9, length.out = ifelse(length(samples) <
                                                                                                                                                          4, 15, 8)), yend = seq(0, 0.9, length.out = ifelse(length(samples) < 4, 15, 8)) + 0.1, x = loc[as.character(.data$GENE)] - 0.49, xend = loc[as.character(.data$GENE)] + 0.49)
      haplo_db_clust_texture <- bind_rows(haplo_db_clust_texture, tmp_point)
    }
  }


  allele_cols <- gsub("IG[H|K|L]", "", gsub("_", "*", hapBy_cols))
  # Adding white space to plot
  heatmap.df <- c()
  samples_order <- c()
  samples_label <- c()
  for (i in 1:length(samples)) {
    samp = samples[i]
    sub <- haplo_db_clust[haplo_db_clust$SUBJECT == samp, ]
    sub2 <- sub[sub$hapBy == allele_cols[1], ]
    sub2$freq <- 0
    sub2$grouper_y <- paste0(sub2$SUBJECT[1], " NA")
    if (i != length(samples)) {
      sub <- rbind(sub, sub2)
      heatmap.df <- rbind(heatmap.df, sub)
      samples_order <- c(samples_order, unique(sub$grouper_y))
      tmp_l <- c(unique(sub$hapBy), "")
      names(tmp_l) <- unique(sub$grouper_y)
      samples_label <- c(samples_label, tmp_l)
    } else {
      heatmap.df <- rbind(heatmap.df, sub)
      samples_order <- c(samples_order, unique(sub$grouper_y))
      tmp_l <- unique(sub$hapBy)
      names(tmp_l) <- unique(sub$grouper_y)
      samples_label <- c(samples_label, tmp_l)
    }
  }

  heatmap.df$grouper_y <- factor(heatmap.df$grouper_y, levels = samples_order)
  if (length(grep("[0-9][0-9]_[0-9][0-9]", heatmap.df$ALLELES)) != 0) {
    non_reliable_alleles_text <- nonReliableAllelesText(heatmap.df[!grepl("NA", heatmap.df$grouper_y) & grepl("[0-9][0-9]_[0-9][0-9]", heatmap.df$ALLELES),
                                                                   ],size = 3)
    non_reliable_alleles_text$grouper_y <- factor(non_reliable_alleles_text$grouper_y, levels = samples_order)
    heatmap.df$ALLELES[grep("[0-9][0-9]_[0-9][0-9]", heatmap.df$ALLELES)] <- "NRA"
    allele_palette <- alleleHapPalette(heatmap.df$ALLELES)
    non_reliable_alleles_text$ALLELES <- factor(non_reliable_alleles_text$ALLELES, levels = allele_palette$AlleleCol)
  } else {
    non_reliable_alleles_text <- c()
    allele_palette <- alleleHapPalette(heatmap.df$ALLELES, NRA = FALSE)
  }

  heatmap.df$ALLELES <- factor(heatmap.df$ALLELES, levels = allele_palette$AlleleCol)

  heatmap.df$GENE_LOC <-  sapply(heatmap.df$GENE,function(i) {loc[as.character(i)]})

  hap_plot <- ggplot() + geom_col(data = heatmap.df, mapping = aes_string(x = "GENE_LOC", y = "freq", fill = "ALLELES"), position = "fill", width = 0.95,
                                  na.rm = T) + scale_fill_manual(values = alpha(names(allele_palette$AlleleCol), allele_palette$transper), name = "Alleles", drop = FALSE) + facet_grid(grouper_y ~
                                                                                                                                                                                          grouper_x, as.table = FALSE, switch = "y", labeller = labeller(grouper_y = samples_label), drop = FALSE) + scale_y_continuous(expand = c(0, 0)) +# scale_x_discrete(expand = c(0,0))
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, colour = "black"), strip.text.x = element_blank(), strip.text.y = element_text(angle = 180,
                                                                                                                                                                   size = 14), panel.grid = element_blank(), strip.placement = "outside", axis.ticks.y = element_line(colour = "white"), axis.line.y.left = element_blank(),
          axis.text.y = element_blank(), strip.background.y = element_blank(), strip.background.x = element_blank(), panel.spacing.y = unit(0.9, "pt"), legend.position = "bottom",
          axis.title.x = element_text(size = 15, colour = "black"), legend.justification = "center") + labs(y = "", x = "Gene") + guides(fill = guide_legend(nrow = round(length(allele_palette$AlleleCol)/9),
                                                                                                                                                             order = 1, override.aes = list(color = "#DCDCDC"))) + scale_x_continuous(expand = c(0,0),breaks = 1:length(unique(heatmap.df$GENE)[order(match(unique(heatmap.df$GENE), levels(heatmap.df$GENE)))]),labels = unique(heatmap.df$GENE)[order(match(unique(heatmap.df$GENE), levels(heatmap.df$GENE)))], sec.axis = dup_axis(name = ""))


  if (mark_low_lk & nrow(haplo_db_clust_texture) != 0) {
    haplo_db_clust_texture$grouper_y <- factor(haplo_db_clust_texture$grouper_y, levels = samples_order)
    haplo_db_clust_texture <- haplo_db_clust_texture[!duplicated(haplo_db_clust_texture[, c("GENE", "ALLELES", "K", "points", "SUBJECT")]), ]
    haplo_db_clust_texture$col <- "<1"
    # Get Allele legend
    gt1 = ggplotGrob(hap_plot)

    haplo_db_clust_texture$ALLELES <- factor(haplo_db_clust_texture$ALLELES, levels = allele_palette$AlleleCol)

    hap_plot <- hap_plot + geom_segment(data = haplo_db_clust_texture, mapping = aes_string(x = "x", xend = "xend", y = "points", yend = "yend", color = "col"), colour = "white")  +
      guides(color = "none", fill = "none")

    # Get lK legend
    gt2 = getDigLegend(unique(haplo_db_clust_texture$col))

    leg1 = gtable::gtable_filter(gt1, "guide-box")
    leg2 = gtable::gtable_filter(gt2, "guide-box")
    # Combine the legends
    leg <- cbind(leg1[["grobs"]][[1]], leg2[["grobs"]][[1]], size = "first")
    # Insert legend into g1 (or g2)
    gt1$grobs[gt1$layout$name == "guide-box"][[1]] <- leg
    gt1$grobs[gt1$layout$name == "guide-box"][[1]]$layout[3, c("t", "b")] <- gt1$grobs[gt1$layout$name == "guide-box"][[1]]$layout[1, c("t", "b")]
    gt1$grobs[gt1$layout$name == "guide-box"][[1]]$layout[3, c("l", "r")] <- gt1$grobs[gt1$layout$name == "guide-box"][[1]]$layout[1, c("l", "r")] +
      2
    legend <- get_legend(gt1)
  } else legend <- get_legend(hap_plot)

  if (is.data.frame(non_reliable_alleles_text)) {
    non_reliable_alleles_text$GENE_LOC <-  sapply(non_reliable_alleles_text$GENE,function(i) {loc[as.character(i)]})
    hap_plot <- hap_plot + geom_text(data = non_reliable_alleles_text, aes_string(label = "text", x = "GENE_LOC", y = "pos"), angle = 90, size = non_reliable_alleles_text$size)
  }

  dend_hap <- plot_grid(hap_plot + theme(legend.position = "none", plot.background = element_rect(fill = "transparent", colour = NA),
                                         panel.background = element_rect(fill = "transparent", colour = NA), axis.line.x = element_line(color="black")),
                        plt_dendr, nrow = 1, align = "h", axis = "bt", rel_widths = c(1, 0.25))
  plot(plot_grid(dend_hap, legend, ncol = 1, rel_heights = c(1, 0.1)))
}

########################################################################################################
#' Graphical output of double chromosome deletions
#'
#' The \code{plotDeletionsByBinom} function generates a graphical output of the double chromosome deletions in multiple samples.
#'
#'
#' @param    GENE.usage.df        double chromosome deletion summary table. See details.
#' @param    chain                the IG chain: IGH,IGK,IGL. Default is IGH.
#' @param    genes.low.cer        a vector of IGH genes known to be with low certantiny in the binomial test. Default is IGHV3-43 and IGHV3-20
#' @param    genes.dup            a vector of IGH genes known to have a duplicated gene. Default is IGHD4-11 that his duplicate is IGHD4-4 and IGHD5-18 that his duplicate is IGHD5-5
#'
#' @return
#'
#' A double chromosome deletion visualization.
#'
#' @details
#'
#' A \code{data.frame} created by \code{binom_test_deletion}.
#'
#' @examples
#'
#' # Load example data and germlines
#' data(samples_db)
#'
#' # Infering haplotype
#' deletions_db = deletionsByBinom(samples_db);
#' plotDeletionsByBinom(deletions_db)
#'
#' @export
plotDeletionsByBinom <- function(GENE.usage.df, chain = c("IGH", "IGK", "IGL"), genes.low.cer = c("IGHV3-43", "IGHV3-20"), genes.dup = c("IGHD4-11", "IGHD5-18")) {

    if (missing(chain)) {
        chain = "IGH"
    }
    chain <- match.arg(chain)

    if (!("SUBJECT" %in% names(GENE.usage.df))) {
        GENE.usage.df$SUBJECT <- rep("S1", nrow(GENE.usage.df))
    }

    genes_hap <- unique(substr(GENE.usage.df$GENE, 4, 4))
    GENE.loc.tmp <- GENE.loc[[chain]][GENE.loc[[chain]] %in% GENE.usage.df$GENE]

    GENE.usage.df$GENE2 <- factor(gsub(chain, "", GENE.usage.df$GENE), levels = gsub(chain, "", GENE.loc.tmp))

    colvec <- ifelse(GENE.loc.tmp %in% genes.low.cer, "red", ifelse(GENE.loc.tmp %in% genes.dup, "purple", "black"))

    ### gene usage with deletions in population according to binom test
    p.del <- ggplot(GENE.usage.df %>% filter(.data$DELETION != "Non reliable"), aes_string(x = "GENE2", y = "FRAC")) + geom_boxplot(outlier.colour = NA) + geom_jitter(aes_string(x = "GENE2",
        color = "DELETION"), width = 0.25, size = 0.5) + theme(axis.text.y = element_text(size = 16), axis.title = element_text(size = 16), axis.text.x = element_text(size = 14,
        angle = 90, hjust = 1, vjust = 0.5, color = colvec), legend.text = element_text(size = 16), legend.position = "none") + ylab("Fraction") + xlab("") +
        scale_color_manual(name = "", labels = c("Deletion", "No Deletion", "NA"), values = c("blue", "black", "grey40"), drop = T) + guides(color = guide_legend(override.aes = list(size = 5)))

    ### heat map of deletions in population according to binom test
    GENE.usage.df$DELETION <- factor(GENE.usage.df$DELETION, levels = levels(GENE.usage.df$DELETION))
    if (length(levels(GENE.usage.df$DELETION)) < 4)
        lab = c("Deletion", "No Deletion", "NA") else lab = c("Deletion", "No Deletion", "NA", "Non reliable")
    if (length(levels(GENE.usage.df$DELETION)) < 4)
        col_val = c("#6d6d6d", "#ffffff", "#dedede") else col_val = c("#6d6d6d", "#ffffff", "#dedede", "#ffefd5")
    heatmap.plot <- ggplot(data = GENE.usage.df, aes_string(x = "GENE2", y = "SUBJECT")) + geom_tile(aes_string(fill = "DELETION")) + scale_fill_manual(name = "", labels = lab,
        values = col_val, drop = F) + scale_x_discrete(drop = FALSE) + ylab("Subject") + xlab("Gene") + theme(axis.text = element_text(size = 12, colour = "black"), axis.title.y = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, colour = "black"), legend.key = element_rect(colour = "black", size = 0.5, linetype = "solid"),
        legend.text = element_text(size = 14), legend.direction = "horizontal", legend.justification = "center", legend.box.just = "bottom")

    legend <- cowplot::get_legend(heatmap.plot)
    heatmap.plot <- heatmap.plot + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
    p.del <- p.del + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

    comb <- plot_grid(p.del, heatmap.plot, ncol = 1, rel_heights = c(0.15, 0.3), align = "hv")
    plot(plot_grid(comb, legend, nrow = 2, rel_heights = c(1, 0.2)))
}

########################################################################################################
#' Graphical output of single chromosome deletions
#'
#' The \code{deletionHeatmap} function generates a graphical output of the single chromosome deletions in multiple samples.
#'
#'
#' @param    hap_table            haplotype summary table. See details.
#' @param    chain                the IG chain: IGH,IGK,IGL. Default is IGH.
#' @param    kThreshDel           the minimum lK (log10 of the Bayes factor) used in \code{createFullHaplotype} to call a deletion. Indicates the color for strong deletion. Defualt is 3.
#'
#' @return
#'
#' A single chromosome deletion visualization.
#'
#' @details
#'
#' A \code{data.frame} created by \code{createFullHaplotype}.
#'
#' @examples
#' # Plotting single choromosme deletion from haplotype inference
#' deletionHeatmap(samplesHaplotype)
#' @export
# Not in use, html_output = FALSE. @param    html_output          If TRUE, a html5 interactive graph is outputed insteaed of the normal plot. Defualt is FALSE
deletionHeatmap <- function(hap_table, kThreshDel = 3, chain = c("IGH", "IGK", "IGL")) {

  if (missing(chain)) {
    chain = "IGH"
  }
  chain <- match.arg(chain)

  if (!("SUBJECT" %in% names(hap_table))) {
    hap_table$SUBJECT <- rep("S1", nrow(hap_table))
  }

  hapBy_cols = names(hap_table)[c(3,4)]

  genes_hap <- unique(substr(hap_table$GENE, 4, 4))

  GENE.loc.tmp <- GENE.loc[[chain]][grep(paste0(genes_hap, collapse = "|"), GENE.loc[[chain]])]

  GENE.loc.tmp <- GENE.loc.tmp[GENE.loc.tmp %in% unique(hap_table$GENE)]

  ALLELE_01_col = hapBy_cols[1]
  ALLELE_02_col = hapBy_cols[2]

  ALLELE_01_num = strsplit(ALLELE_01_col, "_")[[1]][2]
  ALLELE_02_num = strsplit(ALLELE_02_col, "_")[[1]][2]

  ### create deletion streches heatmap
  hap_table$K1[is.na(hap_table$K1)] <- 0
  hap_table$K2[is.na(hap_table$K2)] <- 0

  hap_table.del.heatmap <- hap_table %>% rowwise %>% mutate(K = max(as.numeric(.data$K1), as.numeric(.data$K2), na.rm = T)) %>% select_(.dots = c("SUBJECT", "GENE",
                                                                                                                                                  ALLELE_01_col, "K"))

  hap_table.del.heatmap$HapBy <- rep(ALLELE_01_num, nrow(hap_table.del.heatmap))
  names(hap_table.del.heatmap)[3] <- ALLELE_02_col

  hap_table.del.heatmap <- rbind(hap_table.del.heatmap, data.frame(hap_table %>% rowwise %>% mutate(K = max(as.numeric(.data$K1), as.numeric(.data$K2), na.rm = T)) %>%
                                                                     select_(.dots = c("SUBJECT", "GENE", ALLELE_02_col, "K")), HapBy = ALLELE_02_num))

  names(hap_table.del.heatmap)[3] <- "ALLELE"

  hap_table.del.heatmap$K[hap_table.del.heatmap$K == Inf] <- 0
  hap_table.del.heatmap$K[hap_table.del.heatmap$K == -Inf] <- 0

  hap_table.del.heatmap$DEL <- ifelse(hap_table.del.heatmap$ALLELE == "Del" & hap_table.del.heatmap$K >= kThreshDel, 3, 0)
  hap_table.del.heatmap$DEL[(!hap_table.del.heatmap$ALLELE %in% c("Del", "Unk", "NA")) & hap_table.del.heatmap$K < 3] <- 1
  hap_table.del.heatmap$DEL[hap_table.del.heatmap$ALLELE == "Del" & hap_table.del.heatmap$K < kThreshDel] <- 2
  hap_table.del.heatmap$DEL[hap_table.del.heatmap$ALLELE == "NA"] <- 4



  # manual reshape
  hap_table.del.heatmap.02 <- hap_table.del.heatmap[hap_table.del.heatmap$HapBy == ALLELE_01_num, ]
  hap_table.del.heatmap.03 <- hap_table.del.heatmap[hap_table.del.heatmap$HapBy == ALLELE_02_num, ]
  hap_table.del.heatmap.02$GENE2 <- gsub("IG[H|K|L]", "", hap_table.del.heatmap.02$GENE)
  hap_table.del.heatmap.03$GENE2 <- gsub("IG[H|K|L]", "", hap_table.del.heatmap.03$GENE)

  hap_table.del.heatmap.02$SUBJECT <- factor(x = hap_table.del.heatmap.02$SUBJECT, levels = unique(hap_table.del.heatmap.02$SUBJECT))

  hap_table.del.heatmap.02$GENE2 <- factor(x = hap_table.del.heatmap.02$GENE2, levels = gsub("IG[H|K|L]", "", GENE.loc.tmp))

  hap_table.del.heatmap.03$SUBJECT <- factor(x = hap_table.del.heatmap.03$SUBJECT, levels = unique(hap_table.del.heatmap.03$SUBJECT))

  hap_table.del.heatmap.03$GENE2 <- factor(x = hap_table.del.heatmap.03$GENE2, levels = gsub("IG[H|K|L]", "", GENE.loc.tmp))
  heatmap.df <- rbind(hap_table.del.heatmap.02, hap_table.del.heatmap.03)

  ALLELE_01_col = gsub("_", "*", gsub("IG[H|K|L]", "", ALLELE_01_col))
  ALLELE_02_col = gsub("_", "*", gsub("IG[H|K|L]", "", ALLELE_02_col))
  heatmap.df$DEL <- factor(heatmap.df$DEL, levels = c(0:4))
  heatmap.df$HapBy <- ifelse(heatmap.df$HapBy == ALLELE_01_num, ALLELE_01_col, ALLELE_02_col)
  heatmap.plot <- ggplot(data = heatmap.df, aes_string(x = "GENE2", y = "SUBJECT")) + theme_bw() + geom_tile(aes_string(fill = "DEL")) + facet_wrap(~HapBy, nrow = 2) + scale_x_discrete(drop = FALSE) +
    scale_fill_manual(name = "lK", labels = c("No deletion (lK>=3)", "No deletion (lK<3)", paste0("Deletion (lK<", kThreshDel, ")"), paste0("Deletion (lK>=",
                                                                                                                                            kThreshDel, ")"), "NA"), values = c("white", "#ffb6c1", "lightblue", "#6d6d6d", "#dedede"), drop = FALSE) + ylab("Subject") + xlab("Gene") + theme(strip.text = element_text(size = 18), strip.background =element_rect(fill="seashell2"),
                                                                                                                                                                                                                                                                                             axis.title = element_text(size = 18), axis.text = element_text(size = 14, colour = "black"), axis.text.x = element_text(angle = 90, vjust = 0.5,
                                                                                                                                                                                                                                                                                                                                                                                                                     hjust = 1), plot.margin = margin(b = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
                                                                                                                                                                                                                                                                                             legend.direction = "horizontal", legend.justification = "center", legend.position = "bottom", legend.text = element_text(size = 14), legend.title = element_text(size = 14),
                                                                                                                                                                                                                                                                                             legend.key = element_rect(fill = "white", colour = "black"), axis.line = element_line(colour = "black"))


  del.df.heatmap <- heatmap.df
  del.df.heatmap <- del.df.heatmap %>% filter(.data$ALLELE == "Del")

  del.df.heatmap.cnt <- del.df.heatmap %>% ungroup() %>% group_by(.data$SUBJECT, .data$GENE2) %>% mutate(n = n()) %>% dplyr::slice(1)
  del.df.heatmap.cnt$HapBy[del.df.heatmap.cnt$n == 2] <- "Both"
  del.df.heatmap.cnt <- del.df.heatmap.cnt %>% ungroup() %>% group_by(.data$GENE2, .data$HapBy) %>% count_()
  names(del.df.heatmap.cnt)[3] <- 'n'

  del.df.heatmap.cnt$HapBy <- factor(del.df.heatmap.cnt$HapBy, levels = c(ALLELE_01_col, ALLELE_02_col, "Both"))
  pdel <- ggplot(del.df.heatmap.cnt, aes_string(x = "GENE2", y = "n", fill = "HapBy")) + theme_bw() + geom_bar(stat = "identity", position = "stack", na.rm = T) +
    theme(strip.background = element_blank(), axis.text = element_text(size = 14, colour = "black"), axis.text.x = element_text(angle = 90, vjust = 0.5,
                                                                                                                                hjust = 1), plot.margin = margin(0, 8, 0, 7, "pt"), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), legend.background = element_blank(), panel.background = element_blank()) + ylab("Number of individuals\nwith a deletion") +
    scale_fill_manual(name = "Chromosome", values = c("darksalmon", "deepskyblue", "darkolivegreen3", "grey50")) + scale_x_discrete(drop = FALSE) + xlab("")

  pdel <- pdel + theme(legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center", panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key = element_rect(fill = "white",
                                                                                                                                                                     colour = "black"), legend.margin = margin(0, 0, 0, 0), legend.text = element_text(size = 14), legend.title = element_text(size = 14))

  #a start for html output
  # if (html_output) {
  #
  #     pdel.l <- ggplotly(pdel, height = 400, width = 1500)
  #
  #     heatmap.plot.l <- ggplotly(heatmap.plot, height = 900, width = 1500)
  #
  #     heatmap.bar.l <- subplot(pdel.l, heatmap.plot.l, nrows = 2, shareX = F, heights = c(0.2, 0.6), margin = 0.08)
  #
  #     # fix header box of facet_grif plots
  #
  #     heatmap.bar.l$x$layout$annotations[[2]]$font$size <- 17
  #     heatmap.bar.l$x$layout$annotations[[3]]$font$size <- 17
  #     heatmap.bar.l$x$layout$shapes[[3]]$yanchor <- heatmap.bar.l$x$layout$annotations[[2]]$y
  #     heatmap.bar.l$x$layout$shapes[[5]]$yanchor <- heatmap.bar.l$x$layout$annotations[[3]]$y
  #
  # }

  plot(plot_grid(pdel, heatmap.plot, ncol = 1, rel_heights = c(0.15, 0.35), align = "hv", axis = "b"))
}

##########################################################################
#' Graphical output for single chromosome D or J gene deletions according to V pooled method
#'
#' The \code{plotDeletionsByVpooled} function generates a graphical output for single chromosome D or J gene deletions (for heavy chain only).
#'
#'
#' @param  del.df   a \code{data.frame} created by \code{deletionsByVpooled}.
#' @param  K_ranges vector of one or two integers for log(K) certainty level thresholds
#'
#' @return
#'
#' A single chromosome deletion visualization.
#'
#' @details
#'
#' A \code{data.frame} created by \code{deletionsByVpooled}.
#'
#' @examples
#' \donttest{
#' # Load example data and germlines
#' data(samples_db)
#' del_db <- deletionsByVpooled(samples_db)
#' plotDeletionsByVpooled(del_db)
#' }
#' @export
plotDeletionsByVpooled <- function(del.df, K_ranges = c(3, 7)) {


    if (!("SUBJECT" %in% names(del.df))) {
        del.df$SUBJECT <- rep("S1", nrow(del.df))
    }

    if (length(K_ranges) == 2) {
        del.df$EVENT <- unlist(sapply(1:nrow(del.df), function(i) {
            if (del.df$DELETION[i] > 0 & del.df$K[i] < K_ranges[1])
                return(1)
            if (del.df$DELETION[i] > 0 & del.df$K[i] >= K_ranges[1] & del.df$K[i] < K_ranges[2])
                return(2)
            if (del.df$DELETION[i] > 0 & del.df$K[i] > K_ranges[2])
                return(3)
            if (del.df$DELETION[i] == 0 & del.df$K[i] < K_ranges[1])
                return(4)
            if (del.df$DELETION[i] == 0 & del.df$K[i] >= K_ranges[1] & del.df$K[i] < K_ranges[2])
                return(5)
            if (del.df$DELETION[i] == 0 & del.df$K[i] > K_ranges[2])
                return(6)
        }))


        del.df$EVENT <- factor(del.df$EVENT, levels = 1:6)
        del.df$GENE2 <- factor(gsub("IGH", "", del.df$GENE), levels = gsub("IGH", "", GENE.loc[["IGH"]]))
        labels1 = c(paste0("Deletion lK<", K_ranges[1]), paste0("Deletion ", K_ranges[1], "<=lK<", K_ranges[2]), paste0("Deletion lK>=", K_ranges[2]), paste0("No deletion lK<",
            K_ranges[1]), paste0("No deletion ", K_ranges[1], "<=lK<", K_ranges[2]), paste0("No deletion lK>=", K_ranges[2]))
        values1 = c("lightblue", "cornflowerblue", "#6d6d6d", "lightpink", "lightcoral", "white")
    }

    if (length(K_ranges) == 1) {
        del.df$EVENT <- unlist(sapply(1:nrow(del.df), function(i) {
            if (del.df$DELETION[i] > 0 & del.df$K[i] < K_ranges[1])
                return(1)
            if (del.df$DELETION[i] > 0 & del.df$K[i] >= K_ranges[1])
                return(2)
            if (del.df$DELETION[i] == 0 & del.df$K[i] < K_ranges[1])
                return(3)
            if (del.df$DELETION[i] == 0 & del.df$K[i] >= K_ranges[1])
                return(4)
        }))


        del.df$EVENT <- factor(del.df$EVENT, levels = 1:4)
        del.df$GENE2 <- factor(gsub("IGH", "", del.df$GENE), levels = gsub("IGH", "", GENE.loc[["IGH"]]))
        labels1 = c(paste0("Deletion lK<", K_ranges[1]), paste0("Deletion lK>=", K_ranges[1]), paste0("No deletion lK<", K_ranges[1]), paste0("No deletion lK>=",
            K_ranges[1]))
        values1 = c("lightblue", "#6d6d6d", "lightpink", "lightgrey")
    }

    heatmap.plot <- ggplot(data = del.df, aes_string(x = "GENE2", y = "SUBJECT")) + geom_tile(aes_string(fill = "EVENT")) + scale_fill_manual(name = "", labels = labels1, values = values1,
        drop = F) + ylab("Subject") + xlab("Gene") + theme(strip.text = element_text(size = 14), axis.title = element_text(size = 14), axis.text = element_text(size = 12, colour = "black"),
        legend.position = "bottom", legend.justification = "center", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.key = element_rect(colour = "black",
            size = 0.5, linetype = "solid"), legend.text = element_text(size = 14), axis.line.x.bottom = element_line(colour = "black", inherit.blank = F),
        axis.line.y.left = element_line(colour = "black", inherit.blank = F), axis.line = element_blank(), panel.grid = element_blank(), panel.background = element_blank())



    plot(heatmap.plot)

}
