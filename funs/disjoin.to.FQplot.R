
disjoin.to.FQplot <- function (list_segs, gain, loss) {

  filt_list <- list_segs

  all_segs_filt <- as.data.frame(do.call(rbind, lapply(filt_list,function(x){x[,1:5]})))
  all_segs_filt$chr <- paste("chr", all_segs_filt$chr, sep="")

  all_segs_filt_by_chr <- split(all_segs_filt, all_segs_filt$chr)


  disJoint_by_chr <- lapply(1:length(all_segs_filt_by_chr), function(x){

    chrom <- names(all_segs_filt_by_chr)[x]
    mat <- all_segs_filt_by_chr[[x]]
    z <- GenomicRanges::GRanges(seqnames=mat$chr,
      ranges=IRanges::IRanges(start=mat$loc.start, end=mat$loc.end),
      "sample"=mat$ID, "segmean"=mat$seg.mean)

      disJoint_z <- GenomicRanges::disjoin(z, with.revmap=TRUE)
      disJoint_df <- as.data.frame(disJoint_z)
      colnames(disJoint_df)[1] <- "chr"

      dd <- as.data.frame(do.call(rbind, lapply(1:nrow(disJoint_df), function(x){
        # row <- disJoint_df[x,1:4]
        iters <- disJoint_df[x,6][[1]]
        events <- as.vector(do.call( cbind,lapply(iters, function(x){
          seg.mean <- as.numeric(as.character(mat[x,5]))
          ans <- "normal"
          if (seg.mean>=gain) {
            ans <- "gain"
          }else if (seg.mean<=loss) {
            ans <- "loss"
          }
        })))
        n_gains <- length(which(events=="gain"))
        n_losses <- length(which(events=="loss"))
        n_normals <- length(which(events=="normal"))
        new_row <- cbind(disJoint_df[x,1:4], n_gains, n_losses, n_normals)
        # new_row[1,1] <- chrom
        new_row
      })))
    })
    names(disJoint_by_chr) <- names(all_segs_filt_by_chr)

    sort.chroms <- function(x){
      seq_chroms <- paste("chr", seq(1:length(x)), sep="")
      # print(seq_chroms)
      ans <- rep(NA, length(seq_chroms))
      for (t in 1:length(x)){
        item <- x[t]
        p <- which(seq_chroms==item)
        ans[p] <- t
      }
      ans
    }
    re_order <- sort.chroms(names(disJoint_by_chr))

    disJoint_by_chr <- disJoint_by_chr[re_order]
    fq_df <- as.data.frame(do.call(rbind, disJoint_by_chr))

    l_samples <- length(unique(all_segs_filt$ID))

    fq_df2 <- data.frame(fq_df[,1:3], round((fq_df[5:7]/l_samples)*100,2))
    colnames(fq_df2)[4:6] <- c("gain", "loss", "normal")
    fq_df2
  }
