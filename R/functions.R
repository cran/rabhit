# RAbHIT functions -----------------------------------------------------

#' @include rabhit.R
#' @include internal_functions.R
NULL

##########################################################################
#' Anchor gene haplotype inference
#'
#' The \code{createFullHaplotype} functions infers haplotype based on an anchor gene.
#'
#'
#' @param    clip_db               a \code{data.frame} in Change-O format. See details.
#' @param    toHap_col             a vector of column names for which a haplotype should be inferred. Default is V_CALL and D_CALL.
#' @param    hapBy_col             column name of the anchor gene. Default is J_CALL.
#' @param    hapBy                 a string of the anchor gene name. Default is IGHJ6.
#' @param    toHap_GERM            a vector of named nucleotide germline sequences matching the allele calls in \code{toHap_col} columns in clip_db.
#' @param    relative_freq_priors  if TRUE, the priors for Bayesian inference are estimated from the relative frequencies in clip_db. Else, priors are set to \code{c(0.5,0.5)}. Defualt is TRUE
#' @param    kThreshDel            the minimum lK (log10 of the Bayes factor) to call a deletion. Defualt is 3.
#' @param    rmPseudo              if TRUE non-functional and pseudo genes are removed. Defualt is TRUE.
#' @param    deleted_genes         double chromosome deletion summary table. A \code{data.frame} created by \code{deletionsByBinom}.
#' @param    nonReliable_Vgenes    a list of known non reliable gene assignments. A \code{list} created by \code{nonReliableVGenes}.
#' @param    min_minor_fraction    the minimum minor allele fraction to be used as an anchor gene. Default is 0.3
#' @param    chain                 the IG chain: IGH,IGK,IGL. Default is IGH.
#'
#' @return
#' A \code{data.frame}, in which each row is the haplotype inference summary of a gene from the column selected in \code{toHap_col}.
#'
#'The output containes the following columns:
#' \itemize{
#'  \item \code{SUBJECT}:        the subject name.
#'  \item \code{GENE}:           the gene name.
#'  \item Anchor gene allele 1:  the haplotype inference for chromosome one. The column name is the anchor gene with the first allele.
#'  \item Anchor gene allele 2:  the haplotype inference for chromosome two. The column name is the anchor gene with the second allele.
#'  \item \code{ALLELES}:        allele calls for the gene.
#'  \item \code{PRIORS_ROW}:     priors based on relative allele usage of the anchor gene.
#'  \item \code{PRIORS_COL}:     priors based on relative allele usage of the inferred gene.
#'  \item \code{COUNTS1}:        the appereance count on each chromosome of the first allele from \code{ALLELES}, the counts are seperated by a comma.
#'  \item \code{K1}:             the Bayesian factor value for the first allele (from \code{ALLELES}) inference.
#'  \item \code{COUNTS2}:        the appereance count on each chromosome of the second allele from \code{ALLELES}, the counts are seperated by a comma.
#'  \item \code{K2}:             the Bayesian factor value for the second allele (from \code{ALLELES}) inference.
#'  \item \code{COUNTS3}:        the appereance count on each chromosome of the third allele from \code{ALLELES}, the counts are seperated by a comma.
#'  \item \code{K3}:             the Bayesian factor value for the third allele (from \code{ALLELES}) inference.
#'  \item \code{COUNTS4}:        the appereance count on each chromosome of the fourth allele from \code{ALLELES}, the counts are seperated by a comma.
#'  \item \code{K4}:             the Bayesian factor value for the fourth allele (from \code{ALLELES}) inference.
#'}
#'
#' @details
#' Function accepts a \code{data.frame} in Change-O format (\url{https://changeo.readthedocs.io/en/version-0.4.1---airr-standards/standard.html}) containing the following columns:
#' \itemize{
#'   \item \code{'SUBJECT'}: The subject name
#'   \item \code{'V_CALL'}: V allele call(s) (in an IMGT format)
#'   \item \code{'D_CALL'}: D allele call(s) (in an IMGT format, only for heavy chains)
#'   \item \code{'J_CALL'}: J allele call(s) (in an IMGT format)
#' }
#'
#' @examples
#' # Load example data and germlines
#' data(samples_db, HVGERM, HDGERM)
#'
#' # Selecting a single individual
#' clip_db = samples_db[samples_db$SUBJECT=='I5', ]
#'
#' # Infering haplotype
#' haplo_db = createFullHaplotype(clip_db,toHap_col=c('V_CALL','D_CALL'),
#' hapBy_col='J_CALL',hapBy='IGHJ6',toHap_GERM=c(HVGERM,HDGERM))
#'
#'
#' @export

createFullHaplotype <- function(clip_db, toHap_col = c("V_CALL", "D_CALL"), hapBy_col = "J_CALL", hapBy = "IGHJ6", toHap_GERM, relative_freq_priors = TRUE,
                                kThreshDel = 3, rmPseudo = TRUE, deleted_genes = c(), nonReliable_Vgenes = c(), min_minor_fraction = 0.3, chain = c("IGH", "IGK", "IGL")) {


  # Check if germline was inputed
  if (missing(toHap_GERM))
    stop("Missing toHap_GERM, please input germline sequences")

  if (missing(chain)) {
    chain = "IGH"
  }
  chain <- match.arg(chain)

  if (!("SUBJECT" %in% names(clip_db))) {

    clip_db$SUBJECT <- "S1"
  }

  haplo_db <- c()
  clip_db <- clip_db %>% select(.data$SUBJECT, !!eval(c(hapBy_col,toHap_col)))
  for (sample_name in unique(clip_db$SUBJECT)) {

    if (is.list(nonReliable_Vgenes)){
      nonReliable_Vgenes_vec <- nonReliable_Vgenes[[sample_name]]
    } else nonReliable_Vgenes_vec <- nonReliable_Vgenes

    if (is.data.frame(deleted_genes)) {
      deleted_genes_vec <- deleted_genes %>% filter(.data$SUBJECT == sample_name, .data$DELETION == "Deletion") %>% select(.data$GENE) %>% pull()
      if (is.null(nonReliable_Vgenes_vec))
        nonReliable_Vgenes_vec <- deleted_genes %>% filter(.data$SUBJECT == sample_name, .data$DELETION == "Non reliable") %>% select(.data$GENE) %>% pull()
    } else deleted_genes_vec <- c()


    ### Check if haplotype can be infered by the specific gene in the data set.  Only relevant genes with one assignment
    clip_db_sub <- clip_db %>% filter(.data$SUBJECT==sample_name, !grepl(',', !!as.name(hapBy_col), perl = T))
    hapBy_priors <- clip_db_sub %>% filter(grepl(paste0(hapBy, "\\*"), !!as.name(hapBy_col), perl = T)) %>% dplyr::count(!!as.name(hapBy_col)) %>% mutate(freq = n / sum(n)) %>% select(-n)
    hapBy_priors <- setNames(hapBy_priors[[2]], hapBy_priors[[1]])
    hapBy_alleles <- names(hapBy_priors)
    if (length(hapBy_alleles) != 2)
      stop("Can not haplotype by more or less than two alleles")
    if (min(hapBy_priors) < min_minor_fraction)
      stop("Can not haplotype, minor allele fraction lower than the cutoff set by the user")

    GENES <- unique(gsub("\\*.*", "\\1", grep("^(?=.*IG)", unique(unlist(clip_db_sub[clip_db_sub[, hapBy_col] %in% hapBy_alleles,toHap_col], use.names = F)), value = T, perl = T), perl = T))

    GENES.ref <- unique(gsub("\\*.*", "\\1", names(toHap_GERM), perl = T))

    if (rmPseudo) {
      GENES <- GENES[!grepl("OR|NL", GENES)]
      GENES <- GENES[!(GENES %in% PSEUDO[[chain]])]

      GENES.ref <- GENES.ref[!grepl("OR|NL", GENES.ref)]
      GENES.ref <- GENES.ref[!(GENES.ref %in% PSEUDO[[chain]])]
    }

    GENES.df.num <- data.table::rbindlist(lapply(intersect(GENES, GENES.ref),function(G){

      if(G %in% deleted_genes_vec || G %in% nonReliable_Vgenes_vec){

        relFreqDf.tmp <- data.frame(matrix(c(sample_name, G,
                                             rep(ifelse(G %in% nonReliable_Vgenes_vec, 'NR',
                                                        ifelse(G %in% deleted_genes_vec, 'Del')), 2),
                                             rep(NA, 11)),nrow = 1),stringsAsFactors = F)

        relFreqDf.tmp <- asNum(relFreqDf.tmp)
        return(relFreqDf.tmp)
      } else{

        toHap_col_tmp <-   toHap_col[stringi::stri_detect_fixed(pattern = substr(G,4,4),str = toHap_col)]

        clip_db_sub.G <- clip_db_sub %>% filter(grepl(paste0("^(",G,"\\*[[:digit:]]*[\\_[[:alnum:]]*]*,?)+$"),!!as.name(toHap_col_tmp),perl = T))
        if(substr(G,4,4)=='V'){
          tmp <- clip_db_sub.G %>% filter(stringi::stri_detect_regex(pattern = ",",str = !!as.name(toHap_col_tmp), negate = T))
          tmp2 <- data.table(clip_db_sub.G %>% filter(stringi::stri_detect_regex(pattern = ",",str = !!as.name(toHap_col_tmp))), key=toHap_col_tmp)
          tmp2 <- tmp2[, n:=.N, by=eval(toHap_col_tmp)][,"prop" := n/nrow(clip_db_sub.G)][get("prop") > 0.4,]
          if(length(tmp2[[eval(toHap_col_tmp)]])>0) tmp2[[toHap_col_tmp]] <- alleleCollapse(tmp2[[eval(toHap_col_tmp)]])
          clip_db_sub.G <- rbind(tmp,as.data.frame(tmp2[, c("n", "prop"):=NULL]))

        } else clip_db_sub.G <-  clip_db_sub.G %>% filter(!grepl(',', !!as.name(toHap_col_tmp)))

        tmp <- clip_db_sub.G %>% filter(grepl(paste0('^(?=.*',hapBy,'\\*)'), !!as.name(hapBy_col), perl = T)) %>% select(!!as.name(toHap_col_tmp), !!as.name(hapBy_col)) %>% table()

        if(nrow(tmp)==0){return()
        }else{
          # if one column add the second
          if (ncol(tmp) == 1) {
            toadd <- setdiff(hapBy_alleles, colnames(tmp))
            tmp <- cbind(tmp, rep(0, nrow(tmp)))
            colnames(tmp)[2] <- toadd

            tmp <- as.data.frame(tmp)
            tmp <- tmp[order(colnames(tmp))]
            tmp <- as.matrix(tmp)
          }

          if (relative_freq_priors) {

            clip_db_sub.hapBy <- clip_db_sub.G[clip_db_sub.G[, toHap_col_tmp] %in% rownames(tmp), ]
            toHap_priors <- table(clip_db_sub.hapBy[, toHap_col_tmp])/sum(table(clip_db_sub.hapBy[, toHap_col_tmp]))

            if (length(toHap_priors) != nrow(tmp)) {
              toHap_priors_tmp <- c(rep(0, nrow(tmp)))
              names(toHap_priors_tmp) <- rownames(tmp)
              for (i in names(toHap_priors)) {
                toHap_priors_tmp[i] <- toHap_priors[i]
              }

              toHap_priors <- toHap_priors_tmp
            }

            hap.df <- createHaplotypeTable(tmp, HapByPriors = hapBy_priors, toHapByCol = TRUE, toHapPriors = toHap_priors)

            relFreqDf.tmp <- data.frame(c(sample_name, G, hap.df[, 2:length(hap.df)]),
                                        stringsAsFactors = F)
            relFreqDf.tmp <- asNum(relFreqDf.tmp)
            return(relFreqDf.tmp)
          } else {
            hap.df <- createHaplotypeTable(tmp)
            relFreqDf.tmp <- data.frame(c(sample_name, G, hap.df[, 2:length(hap.df)]),
                                        stringsAsFactors = F)
            relFreqDf.tmp <- asNum(relFreqDf.tmp)
            return(relFreqDf.tmp)

          }
        }
      }
    }
    ),use.names=FALSE) %>% as.data.frame()

    colnames(GENES.df.num) <- c("SUBJECT", "GENE", gsub(pattern = "*", "_", hapBy_alleles, fixed = T),
                                'ALLELES', 'PRIORS_ROW', 'PRIORS_COL','COUNTS1', 'K1',
                                'COUNTS2', 'K2','COUNTS3', 'K3','COUNTS4', 'K4')
    # Check if toHap_col genes are in toHap_GERM
    if (length(GENES.df.num) == 0)
      stop("Genes in haplotype column to be infered do not match the genes germline given")


    ## Fill deleted according to a k thershold
    unkIDX <- which(GENES.df.num[gsub("*", "_", hapBy_alleles[1], fixed = T)][, 1] == "Unk")
    delIDX <- unkIDX[which(sapply(unkIDX, function(i) min(GENES.df.num[i,paste0('K',1:4)], na.rm = T)) >= kThreshDel)]
    GENES.df.num[delIDX, paste(gsub("*", "_", hapBy_alleles[1], fixed = T))] <- "Del"

    unkIDX <- which(GENES.df.num[gsub("*", "_", hapBy_alleles[2], fixed = T)][, 1] == "Unk")
    delIDX <- unkIDX[which(sapply(unkIDX, function(i) min(GENES.df.num[i,paste0('K',1:4)], na.rm = T)) >= kThreshDel)]
    GENES.df.num[delIDX, paste(gsub("*", "_", hapBy_alleles[2], fixed = T))] <- "Del"

    ## Add as unknown genes that do not appear in the individual and mark them as unknown

    GENES.MISSING <- GENES.ref[!(GENES.ref %in% GENES.df.num$GENE)]
    if (length(GENES.MISSING) > 0) {
      m <- length(GENES.MISSING)

      sub.df <- data.frame(do.call(rbind, lapply(1:m, function(x) {
        c(sample_name, NA, "Unk", "Unk", rep(NA, 11))
      })))
      names(sub.df) <- names(GENES.df.num)

      sub.df$GENE <- GENES.MISSING

      sub.df[,gsub("*", "_", hapBy_alleles, fixed = T)] <- matrix(rep(ifelse(GENES.MISSING %in% nonReliable_Vgenes_vec,'NR',ifelse(GENES.MISSING %in% deleted_genes_vec, 'Del', 'Unk')), 2), ncol=2)

      GENES.df.num <- rbind(GENES.df.num, sub.df)
    }

    haplo_db <- rbind(haplo_db, GENES.df.num)
  }

  return(haplo_db)

}

##########################################################################
#' Double chromosome deletion by relative gene usage
#'
#' The \code{deletionsByBinom} function inferes double chromosome deletion events by relative gene usage.
#'
#'
#' @param    clip_db               a \code{data.frame} in Change-O format. See details.
#' @param    chain                 the IG chain: IGH,IGK,IGL. Default is IGH.
#' @param    nonReliable_Vgenes    a list of known non reliable gene assignments. A \code{list} created by \code{nonReliableVGenes}.
#'
#' @return  A \code{data.frame}, in which each row is the double chomosome deletion inference of a gene.
#'
#'The output containes the following columns:
#' \itemize{
#'  \item \code{SUBJECT}:       the subject name.
#'  \item \code{GENE}:          the gene call
#'  \item \code{FRAC}:          the relative gene usage of the gene
#'  \item \code{CUTOFF}:        the the cutoff of for the binomial test
#'  \item \code{PVAL}:          the p-value of the binomial test
#'  \item \code{DELETION}:      if a double chromosome deletion event of a gene occured.
#'}
#'
#' @details
#' The function accepts a \code{data.frame} in Change-O format (\url{https://changeo.readthedocs.io/en/version-0.4.1---airr-standards/standard.html}) containing the following columns:
#' \itemize{
#'   \item \code{'SUBJECT'}: The subject name
#'   \item \code{'V_CALL'}: V allele call(s) (in an IMGT format)
#'   \item \code{'D_CALL'}: D allele call(s) (in an IMGT format, only for heavy chains)
#'   \item \code{'J_CALL'}: J allele call(s) (in an IMGT format)
#' }
#'
#' @examples
#' # Load example data and germlines
#' data(samples_db)
#'
#' # Selecting a single individual
#' clip_db = samples_db[samples_db$SUBJECT=='I5', ]
#' # Infering haplotype
#' del_binom_df = deletionsByBinom(clip_db)
#' head(del_binom_df)
#'
#' @export

deletionsByBinom <- function(clip_db, chain = c("IGH", "IGK", "IGL"), nonReliable_Vgenes = c()) {
    if (missing(chain)) {
        chain = "IGH"
    }
    chain <- match.arg(chain)

    if (!("SUBJECT" %in% names(clip_db))) {
        clip_db$SUBJECT <- "S1"
    }

    GENE.loc.NoPseudo <- GENE.loc[[chain]][!grepl("[OR|NL]", GENE.loc[[chain]])]
    GENE.loc.NoPseudo <- GENE.loc.NoPseudo[!(GENE.loc.NoPseudo %in% PSEUDO[[chain]]) & ( GENE.loc.NoPseudo %in% Binom.test.gene.cutoff[[chain]]$GENE)]

    GENE.usage <- vector("list", length = length(GENE.loc.NoPseudo))
    names(GENE.usage) <- GENE.loc.NoPseudo
    for (samp in unique(clip_db$SUBJECT)) {
        clip_db_sub <- clip_db[clip_db$SUBJECT == samp, ]
        # V gene distribution
        V_CALLS <- table(sapply(strsplit(clip_db_sub$V_CALL, "*", fixed = T), "[", 1))
        V_CALLS_freq <- V_CALLS/sum(V_CALLS)
        for (v in GENE.loc[[chain]]) {
            if (v %in% names(V_CALLS)) {
                GENE.usage[[v]] <- c(GENE.usage[[v]], V_CALLS_freq[v])
                names(GENE.usage[[v]])[length(GENE.usage[[v]])] <- samp
            } else {
                if (grepl(paste0(chain, "V"), v)) {
                  GENE.usage[[v]] <- c(GENE.usage[[v]], 0)
                  names(GENE.usage[[v]])[length(GENE.usage[[v]])] <- samp
                }
            }
        }
        # D gene distribution
        if (chain == "IGH") {
            D_SINGLE <- grep(pattern = ",", clip_db_sub$D_CALL, invert = T)
            D_CALLS <- table(sapply(strsplit(clip_db_sub$D_CALL[D_SINGLE], "*", fixed = T), "[", 1))
            D_CALLS_freq <- D_CALLS/sum(D_CALLS)
            for (d in GENE.loc[[chain]]) {
                if (d %in% names(D_CALLS)) {
                  GENE.usage[[d]] <- c(GENE.usage[[d]], D_CALLS_freq[d])
                  names(GENE.usage[[d]])[length(GENE.usage[[d]])] <- samp
                } else {
                  if (grepl(paste0(chain, "D"), d)) {
                    GENE.usage[[d]] <- c(GENE.usage[[d]], 0)
                    names(GENE.usage[[d]])[length(GENE.usage[[d]])] <- samp
                  }
                }
            }
        }
        # J gene distribution
        J_CALLS <- table(sapply(strsplit(clip_db_sub$J_CALL, "*", fixed = T), "[", 1))
        J_CALLS_freq <- J_CALLS/sum(J_CALLS)
        for (j in GENE.loc[[chain]]) {
            if (j %in% names(J_CALLS)) {
                GENE.usage[[j]] <- c(GENE.usage[[j]], J_CALLS_freq[j])
                names(GENE.usage[[j]])[length(GENE.usage[[j]])] <- samp
            } else {
                if (grepl(paste0(chain, "J"), j)) {
                  GENE.usage[[j]] <- c(GENE.usage[[j]], 0)
                  names(GENE.usage[[j]])[length(GENE.usage[[j]])] <- samp
                }
            }
        }
    }


    ### Gene usage for violin plot for paper
    SAMPLE.SIZE.V <- sapply(unique(clip_db$SUBJECT), function(x) {
        sub <- clip_db[clip_db$SUBJECT == x, ]
        nrow(sub[!grepl(",", sub$V_CALL), ])
    })
    names(SAMPLE.SIZE.V) <- unique(clip_db$SUBJECT)
    if (chain == "IGH") {
        SAMPLE.SIZE.D <- sapply(unique(clip_db$SUBJECT), function(x) {
            sub <- clip_db[clip_db$SUBJECT == x, ]
            nrow(sub[!grepl(",", sub$D_CALL), ])
        })
        names(SAMPLE.SIZE.D) <- unique(clip_db$SUBJECT)
    }
    SAMPLE.SIZE.J <- sapply(unique(clip_db$SUBJECT), function(x) {
        sub <- clip_db[clip_db$SUBJECT == x, ]
        nrow(sub[!grepl(",", sub$J_CALL), ])
    })
    names(SAMPLE.SIZE.J) <- unique(clip_db$SUBJECT)
    SAMPLE.SIZE <- sapply(unique(clip_db$SUBJECT), function(x) nrow(clip_db[clip_db$SUBJECT == x, ]))
    names(SAMPLE.SIZE) <- unique(clip_db$SUBJECT)

    gusage <- unlist(GENE.usage)
    gusage.gene <- sapply(strsplit(names(unlist(GENE.usage)), ".", fixed = T), "[", 1)
    gusage.samp <- sapply(strsplit(names(unlist(GENE.usage)), ".", fixed = T), "[", 2)
    GENE.usage.df <- data.frame(SUBJECT = gusage.samp, GENE = gusage.gene, FRAC = gusage, stringsAsFactors = F, row.names = NULL)

    GENE.usage.df <- GENE.usage.df %>% filter(.data$GENE %in% GENE.loc.NoPseudo)
    GENE.usage.df$NREADS <- SAMPLE.SIZE[GENE.usage.df$SUBJECT]

    GENE.usage.df$min_frac <- sapply(1:nrow(GENE.usage.df), function(x) {
        unique(Binom.test.gene.cutoff[[chain]]$min_frac[Binom.test.gene.cutoff[[chain]]$GENE == GENE.usage.df$GENE[x]])
    })

    GENE.usage.df.V <- GENE.usage.df %>% filter(grepl(paste0(chain, "V"), .data$GENE))
    GENE.usage.df.J <- GENE.usage.df %>% filter(grepl(paste0(chain, "J"), .data$GENE))

    GENE.usage.df.V$NREADS <- SAMPLE.SIZE.V[GENE.usage.df.V$SUBJECT]
    GENE.usage.df.J$NREADS <- SAMPLE.SIZE.J[GENE.usage.df.J$SUBJECT]

    GENE.usage.df.V$NREADS_SAMP <- SAMPLE.SIZE[GENE.usage.df.V$SUBJECT]
    GENE.usage.df.J$NREADS_SAMP <- SAMPLE.SIZE[GENE.usage.df.J$SUBJECT]

    GENE.loc.V <- GENE.loc[[chain]][grep("V", GENE.loc[[chain]])]
    GENE.usage.df.V <- binomTestDeletion(GENE.usage.df.V, cutoff = 0.001, p.val.cutoff = 0.01, chain = chain, GENE.loc.V)

    if (is.list(nonReliable_Vgenes)) {
        for (sample_name in names(nonReliable_Vgenes)) {
            levels(GENE.usage.df.V$col) <- c(levels(GENE.usage.df.V$col), "Non reliable")
            idx <- which(GENE.usage.df.V$GENE[GENE.usage.df.V$SUBJECT == sample_name] %in% nonReliable_Vgenes[[sample_name]])
            GENE.usage.df.V$col[GENE.usage.df.V$SUBJECT == sample_name][idx] <- "Non reliable"
        }
    }

    GENE.loc.J <- GENE.loc[[chain]][grep("J", GENE.loc[[chain]])]
    GENE.usage.df.J <- binomTestDeletion(GENE.usage.df.J, cutoff = 0.005, p.val.cutoff = 0.01, chain = chain, GENE.loc.J)

    if (chain == "IGH") {
        GENE.usage.df.D <- GENE.usage.df %>% filter(grepl(paste0(chain, "D"), .data$GENE))
        GENE.usage.df.D$NREADS <- SAMPLE.SIZE.D[GENE.usage.df.D$SUBJECT]
        GENE.usage.df.D$NREADS_SAMP <- SAMPLE.SIZE[GENE.usage.df.D$SUBJECT]
        GENE.loc.D <- GENE.loc[[chain]][grep("D", GENE.loc[[chain]])]
        GENE.usage.df.D <- binomTestDeletion(GENE.usage.df.D, cutoff = 0.005, p.val.cutoff = 0.01, chain, GENE.loc.D)

        GENE.usage.df.V$col <- factor(GENE.usage.df.V$col, levels = levels(GENE.usage.df.V$col))
        GENE.usage.df.D$col <- factor(GENE.usage.df.D$col, levels = levels(GENE.usage.df.V$col))
        GENE.usage.df.J$col <- factor(GENE.usage.df.J$col, levels = levels(GENE.usage.df.V$col))

        GENE.usage.df.V$GENE <- factor(GENE.usage.df.V$GENE, levels = GENE.loc[[chain]])
        GENE.usage.df.D$GENE <- factor(GENE.usage.df.D$GENE, levels = GENE.loc[[chain]])
        GENE.usage.df.J$GENE <- factor(GENE.usage.df.J$GENE, levels = GENE.loc[[chain]])
        GENE.usage.df <- rbind(GENE.usage.df.V, GENE.usage.df.D, GENE.usage.df.J)
    } else {
        GENE.usage.df.V$col <- factor(GENE.usage.df.V$col, levels = levels(GENE.usage.df.V$col))
        GENE.usage.df.J$col <- factor(GENE.usage.df.J$col, levels = levels(GENE.usage.df.V$col))

        GENE.usage.df.V$GENE <- factor(GENE.usage.df.V$GENE, levels = GENE.loc[[chain]])
        GENE.usage.df.J$GENE <- factor(GENE.usage.df.J$GENE, levels = GENE.loc[[chain]])
        GENE.usage.df <- rbind(GENE.usage.df.V, GENE.usage.df.J)
    }



    GENE.usage.df <- GENE.usage.df %>% ungroup() %>% select(.data$SUBJECT, .data$GENE, .data$FRAC, CUTOFF = .data$min_frac, PVAL = .data$pval_adj, DELETION = .data$col)
    return(GENE.usage.df)
}

##########################################################################
#' Single chromosomal D or J gene deletions inferred by the V pooled method
#'
#' The \code{deletionsByVpooled} function inferes single chromosomal deletion for D and J gene .
#'
#'
#' @param  clip_db                  a \code{data.frame} in Change-O format. See details.
#' @param  chain                    the IG chain: IGH,IGK,IGL. Default is IGH.
#' @param  deletion_col             a vector of column names for which single chromosome deletions should be inferred. Default is J_CALL and D_CALL.
#' @param  count_thresh             integer, the minimun number of sequences mapped to a specific V gene to be included in the V pooled inference.
#' @param  deleted_genes            double chromosome deletion summary table. A \code{data.frame} created by \code{deletionsByBinom}.
#' @param  min_minor_fraction       the minimum minor allele fraction to be used as an anchor gene. Default is 0.3
#' @param  kThreshDel               the minimum lK (log10 of the Bayes factor) to call a deletion. Defualt is 3.
#' @param  nonReliable_Vgenes       a list of known non reliable gene assignments. A \code{list} created by \code{nonReliableVGenes}.
#'
#' @return
#' A \code{data.frame}, in which each row is the single chomosome deletion inference of a gene.
#'
#'The output containes the following columns:
#' \itemize{
#'  \item \code{SUBJECT}:       the subject name.
#'  \item \code{GENE}:          the gene call
#'  \item \code{DELETION}:      chromosome deletions inferred. Encoded 1 for deletion and 0 for no deletion.
#'  \item \code{K}:             the Bayesian factor value for the deletion inference.
#'  \item \code{COUNTS}:        the appereance count of the gene on each chromosome, the counts are seperated by a comma.
#'}
#'
#' @details
#' The function accepts a \code{data.frame} in Change-O format (\url{https://changeo.readthedocs.io/en/version-0.4.1---airr-standards/standard.html}) containing the following columns:
#' \itemize{
#'   \item \code{'SUBJECT'}: The subject name
#'   \item \code{'V_CALL'}: V allele call(s) (in an IMGT format)
#'   \item \code{'D_CALL'}: D allele call(s) (in an IMGT format, only for heavy chains)
#'   \item \code{'J_CALL'}: J allele call(s) (in an IMGT format)
#' }
#'
#' @examples
#' \donttest{
#' data(samples_db)
#'
#' # Infering V pooled deletions
#' del_db <- deletionsByVpooled(samples_db)
#' head(del_db)
#' }
#' @export
# not for light chain
deletionsByVpooled <- function(clip_db, chain = c("IGH","IGK","IGL"), deletion_col = c("D_CALL","J_CALL"), count_thresh = 50, deleted_genes = "", min_minor_fraction = 0.3,
    kThreshDel = 3, nonReliable_Vgenes = c()) {

    if (missing(chain)) {
      chain = "IGH"
    }
    chain <- match.arg(chain)

    if(chain != "IGH") deletion_col = "J_CALL"
    deletion_col <- match.arg(deletion_col)

    if (!("SUBJECT" %in% names(clip_db))) {
        clip_db$SUBJECT <- "S1"
    }

    del.df <- c()
    for (sample_name in unique(clip_db$SUBJECT)) {

        clip_db_sub <- clip_db[clip_db$SUBJECT == sample_name, ]
        if (is.data.frame(deleted_genes))
            deleted_genes_df <- deleted_genes %>% filter(.data$SUBJECT == sample_name,
                                                         grepl(paste0(chain,c("D","J"),collapse = "|"), .data$GENE)) else deleted_genes_df <- c()

        ### remove non reliable V genes prior to pooling
        if (is.list(nonReliable_Vgenes))
            nonReliable_Vgenes_vec <- nonReliable_Vgenes[[sample_name]] else nonReliable_Vgenes_vec <- nonReliable_Vgenes

        ### Only relevant genes with one assignment
        IND <- apply(clip_db_sub, 1, function(x) {
            sum(grep(",", x[c("V_CALL", deletion_col)], invert = F))
        })

        clip_db_sub <- clip_db_sub[IND == 0, ]

        ### Test for heterozygous V genes

        VGENES <- unique(sapply(strsplit(clip_db_sub$V_CALL, split = "*", fixed = T), "[", 1))
        VGENES <- VGENES[!VGENES %in% nonReliable_Vgenes_vec]
        GENES <- unlist(sapply(VGENES, function(x) {
            gene_counts <- table(grep(clip_db_sub$V_CALL, pattern = paste0(x, "*"), fixed = T, value = T))
            if (length(gene_counts) == 2 & (min(gene_counts)/sum(gene_counts)) >= min_minor_fraction & sum(gene_counts) >= count_thresh) {
                return(x)
            }
        }))

        V.df <- list()
        if (length(GENES) > 0) {
            print(paste0("The following genes used for pooled deletion detection for sample ", sample_name))
            print(paste(GENES, sep = ","))
            toHapGerm <- unlist(unname(GERM[[chain]][deletion_col]))
            for (G in GENES) {

                full.hap <- createFullHaplotype(clip_db_sub, toHap_col = deletion_col, hapBy_col = "V_CALL", hapBy = G,
                                                toHap_GERM = toHapGerm, relative_freq_priors = T,
                  kThreshDel = kThreshDel, rmPseudo = T, deleted_genes = deleted_genes_df, chain = chain)

                full.hap$V_ALLELE_1 <- strsplit(names(full.hap)[3], "_")[[1]][2]
                full.hap$V_ALLELE_2 <- strsplit(names(full.hap)[4], "_")[[1]][2]

                names(full.hap)[3] <- paste0(G, "_", 1)
                names(full.hap)[4] <- paste0(G, "_", 2)
                V.df[[G]] <- rbind(V.df[[G]], full.hap)

            }
        } else {
            if(sample_name == tail(unique(clip_db$SUBJECT),1)){
            stop("No heterozygous V genes found for deletion detection, try changing the parameters")
            }else{
              message(paste0("No heterozygous V genes found for ",sample_name," deletion detection, try changing the parameters"))
              next()
            }
        }

        GENES <- unlist(sapply(names(V.df), function(G) {
            if (sample_name %in% V.df[[G]]$SUBJECT) {
                return(G)
            }
        }))

        if (!is.null(GENES)) {
            d.del.df <- c()

            for (G in GENES) {
                tmp <- V.df[[G]] %>% filter(.data$SUBJECT == sample_name, grepl(paste0(chain,c("D","J"),collapse = "|"), .data$GENE))
                tmp$DELETION <- apply(tmp, 1, function(x) {
                  if (x[3] == "Unk" & x[4] != "Unk" | x[3] != "Unk" & x[4] == "Unk") {
                    return(1)
                  }
                  if (x[3] == "Del" & x[4] != "Del" | x[3] != "Del" & x[4] == "Del") {
                    return(2)
                  }
                  if (x[3] == "Unk" & x[4] == "Unk") {
                    return(3)
                  }
                  if (x[3] == "Del" & x[4] == "Del") {
                    return(4)
                  }
                  return(0)
                })
                tmp$K <- sapply(1:nrow(tmp), function(i) {
                  if (!is.na(tmp$COUNTS2[i])) {
                    return(max(tmp$K1[i], tmp$K2[i], na.rm = T))
                  } else {
                    return(tmp$K1[i])
                  }
                })
                tmp$COUNTS <- sapply(1:nrow(tmp), function(i) {
                  if (!is.na(tmp$COUNTS2[i])) {
                    cnt1 <- as.numeric(strsplit(as.character(tmp$COUNTS1[i]), ",")[[1]])
                    cnt2 <- as.numeric(strsplit(as.character(tmp$COUNTS2[i]), ",")[[1]])
                    cnt <- cnt1 + cnt2
                    return(paste0(cnt, collapse = ","))
                  } else {
                    return(as.character(tmp$COUNTS1[i]))
                  }
                })

                tmp <- tmp %>% mutate(SUBJECT = sample_name) %>% select(.data$SUBJECT, .data$GENE, .data$DELETION, .data$K, .data$COUNTS)
                tmp$V_GENE <- rep(G, nrow(tmp))
                tmp$DELETION2 <- ifelse(tmp$DELETION == 0, 0, 1)
                d.del.df <- rbind(d.del.df, tmp)

            }


            tmp.df <- do.call("rbind", lapply(unique(d.del.df$GENE), function(x) {
                x.df <- d.del.df %>% filter(.data$GENE == x) %>% mutate(minCOUNT = (sapply(strsplit(as.character(.data$COUNTS), ","), function(x) {
                  min(as.numeric(x))
                })), maxCOUNT = (sapply(strsplit(as.character(.data$COUNTS), ","), function(x) {
                  max(as.numeric(x))
                })))
                x.df <- x.df %>% filter(!is.na(.data$K))
                if (nrow(x.df) != 0) {

                  x.df <- x.df %>% group_by(.data$DELETION2) %>% mutate(minCOUNTsum = sum(.data$minCOUNT), maxCOUNTsum = sum(.data$maxCOUNT)) %>% slice(1) %>% mutate(DELETION3 = which.max(get_probabilites_with_priors(c(.data$maxCOUNTsum, .data$minCOUNTsum))[1:2]),
                    K2 = max(get_probabilites_with_priors(c(.data$maxCOUNTsum, .data$minCOUNTsum))[1:2]) - min(get_probabilites_with_priors(c(.data$maxCOUNTsum, .data$minCOUNTsum))[1:2]))
                  return(x.df)
                }

            }))


            tmp.df.slct <- tmp.df %>% mutate(COUNTS2 = paste0(.data$minCOUNTsum, ",", .data$maxCOUNTsum)) %>% select(.data$SUBJECT, .data$GENE, .data$DELETION2, .data$DELETION3, .data$K2, .data$COUNTS2, .data$V_GENE)
            tmp.df.slct$V_GENE <- ifelse(tmp.df.slct$DELETION2 != 0, "POOLED_UNK", "POOLED_KNOWN")
            tmp.df.slct <- tmp.df.slct[names(tmp.df.slct) != "DELETION2"]
            d.del.df <- d.del.df[names(d.del.df) != "DELETION2"]
            names(tmp.df.slct) <- names(d.del.df)
            tmp.df.slct$DELETION <- ifelse(tmp.df.slct$DELETION == 2, 0, 1)
            ## ALL Vs
            tmp.df.slct.all <- tmp.df %>% ungroup() %>% group_by(.data$GENE) %>% mutate(minCOUNTsum = sum(.data$minCOUNTsum), maxCOUNTsum = sum(.data$maxCOUNTsum)) %>% slice(1) %>%
                mutate(DELETION3 = which.max(get_probabilites_with_priors(c(.data$maxCOUNTsum, .data$minCOUNTsum))[1:2]), K2 = max(get_probabilites_with_priors(c(.data$maxCOUNTsum,
                .data$minCOUNTsum))[1:2]) - min(get_probabilites_with_priors(c(.data$maxCOUNTsum, .data$minCOUNTsum))[1:2]))
            tmp.df.slct.all <- tmp.df.slct.all %>% mutate(SUBJECT = sample_name, COUNTS = paste0(.data$minCOUNTsum, ",", .data$maxCOUNTsum)) %>% select(.data$SUBJECT, .data$GENE, .data$DELETION3, .data$K2, .data$COUNTS) %>% rename(DELETION = .data$DELETION3, K = .data$K2)
            tmp.df.slct.all$V_GENE <- "V(pooled)"
            tmp.df.slct.all$DELETION <- ifelse(tmp.df.slct.all$DELETION == 2, 0, 1)
            d.del.df <- rbind(d.del.df, tmp.df.slct)
            d.del.df <- rbind(d.del.df, as.data.frame(tmp.df.slct.all))
            d.del.df.pooled <- d.del.df %>% filter(.data$V_GENE == "V(pooled)") %>% select(- .data$V_GENE)
            d.del.df.pooled$K <- round(as.numeric(d.del.df.pooled$K), digits = 2)
            d.del.df.pooled$GENE <- factor(x = d.del.df.pooled$GENE, levels = GENE.loc[[chain]])
        }
        del.df <- rbind(del.df, d.del.df.pooled)
    }
    return(del.df)

}

##########################################################################
#' Detect non reliable gene assignment
#'
#' \code{nonReliableVGenes} Takes a \code{data.frame} in Change-O format  and detect non reliable IGHV genes. A non reliable gene is
#' when the ratio of the multiple assignments with a gene is below the threshold.
#'
#' @param  clip_db              a \code{data.frame} in Change-O format. See details.
#' @param  thresh               the threshold to consider non reliable gene. Defualt is 0.9
#' @param  appearance           the minimun fraction of gene appearance to be considered for reliability check. Defualt is 0.01.
#'
#' @return  a nested list of non reliable genes for all subject.
#'
#' @details
#'
#' The function accepts a \code{data.frame} in Change-O format (\url{https://changeo.readthedocs.io/en/version-0.4.1---airr-standards/standard.html}) containing the following columns:
#' \itemize{
#'   \item \code{'SUBJECT'}: subject names
#'   \item \code{'V_CALL'}: V allele call(s) (in an IMGT format)
#' }
#'
#' @examples
#' # Example IGHV call data frame
#' clip_db <- data.frame(SUBJECT=rep('S1',6),
#' V_CALL=c('IGHV1-69*01','IGHV1-69*01','IGHV1-69*01,IGHV1-69*02',
#' 'IGHV4-59*01,IGHV4-61*01','IGHV4-59*01,IGHV4-31*02','IGHV4-59*01'))
#' # Detect non reliable genes
#' nonReliableVGenes(clip_db)
#' @export
# only heavy chain
nonReliableVGenes <- function(clip_db, thresh = 0.9, appearance = 0.01) {


    # if(missing(chain)) { chain='IGH' } chain <- match.arg(chain)

    if (!("SUBJECT" %in% names(clip_db))) {
        clip_db$SUBJECT <- "S1"
    }
    chain = "IGH"
    GENE.loc.NoPseudo <- GENE.loc[[chain]][grepl("V((?!NL).)*$", GENE.loc[[chain]], perl = T)]
    GENE.loc.NoPseudo <- GENE.loc.NoPseudo[!(GENE.loc.NoPseudo %in% PSEUDO[[chain]])]

    non_reliable_genes <- sapply(as.character(unique(clip_db$SUBJECT)), function(sample_name){
      gene_call <- clip_db$V_CALL[clip_db$SUBJECT == sample_name]
      unlist(lapply(GENE.loc.NoPseudo, function(gene){
        sub <- grep(paste0("((?:,|^)", gene, "\\*[[:digit:]]*[\\_[[:alnum:]]*]*(?:,|$))"), gene_call,perl = T, value = T)
        total <- length(sub)
        if (total/length(gene_call) >= appearance){
          sa <- length(grep(paste0("^(",gene,"\\*[[:digit:]]*[\\_[[:alnum:]]*]*,?)+$"), sub, perl = T))
          if (sa/total < thresh) return(gene)
        }

      }), FALSE, FALSE)}, USE.NAMES = T, simplify = F)
    return(non_reliable_genes)
}
