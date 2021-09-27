seq.pair <- rep(list(c(0, 1)), 6)
seq.pair <- expand.grid(seq.pair)

seq.pair <- seq.pair[rowSums(seq.pair) == 3, ]
rle.seq.pair <- apply(seq.pair, 1, rle)

final.seq.pair <- NULL
for (i in 1:nrow(seq.pair)) {
  
  if (length(rle.seq.pair[[i]]$lengths[rle.seq.pair[[i]]$values == 0]) == 1) {
    next
  }
  
  final.seq.pair <- rbind(final.seq.pair,
                          seq.pair[i, ])
}
colnames(final.seq.pair) <- paste0("p", 1:6)

write.csv(final.seq.pair, file = "TestDesgin.csv", row.names = F)
