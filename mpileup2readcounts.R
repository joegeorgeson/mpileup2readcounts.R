insert.text <- function(string) {
  pattern <- "(?:\\+)(\\d+)(?=[A-Za-z]+)([A-Za-z]+)"
  matches <- gregexpr(pattern, string, perl = TRUE)
  
  result <- lapply(regmatches(string, matches), function(x) {
    num.chars <- unlist(regmatches(x, gregexpr("[0-9]+", x)))
    start.index <- 1+nchar(as.numeric(num.chars))
    end.index <- start.index+as.numeric(num.chars)
    chars <- substr(x,start.index+1, end.index)
    if(length(chars)>0){unlist(chars)}
    else{NA}
  })
  
  unlist(result)
}
del.text <- function(string) {
  pattern <- "(?:\\-)(\\d+)(?=[A-Za-z]+)([A-Za-z]+)"
  matches <- gregexpr(pattern, string, perl = TRUE)
  
  result <- lapply(regmatches(string, matches), function(x) {
    num.chars <- unlist(regmatches(x, gregexpr("[0-9]+", x)))
    start.index <- 1+nchar(as.numeric(num.chars))
    end.index <- start.index+as.numeric(num.chars)
    chars <- substr(x,start.index+1, end.index)
    if(length(chars)>0){unlist(chars)}
    else{NA}
  })
  
  unlist(result)

}
del.ref <- function(string, ref.seq){
  ref.dels <- strsplit(string, "")[[1]][strsplit(string, "")[[1]]=="*"]
  if(length(ref.dels)>0){
    result <- rep(ref.seq, length(ref.dels))
    unlist(result)
  }else{NA}
}
noIndel.text <- function(string) {
  pattern <- "(?:\\-|\\+)(\\d+)(?=[A-Za-z]+)([A-Za-z]+)"
  matches <- gregexpr(pattern, string, perl = TRUE)
  
  to.remove <- regmatches(string, matches)[[1]]
  to.remove.fix <- lapply(to.remove, function(x) {
    num.chars <- unlist(regmatches(x, gregexpr("[0-9]+", x)))
    start.index <- 1+nchar(as.numeric(num.chars))
    end.index <- start.index+as.numeric(num.chars)
    out.fix <- substr(x,start.index+1, start.index+as.numeric(num.chars))
    unlist(out.fix)
  })
  
  string.out <- string
  if(length(to.remove.fix)>0){
    for(i in 1:length(to.remove)){
      string.out <- gsub(to.remove[i], "", string.out)
    }
  }else{string.out <- string}
  return(string.out)
}

table.header <- paste("chr", "coord", "refSeq", "depth", "A", "T", "C", "G", "N",
                  "Starts", "Stops", "Ins", "Del.1bp.DS", "Del.this.pos", sep="\t")

cat(table.header, sep = "\n")

f <- file("stdin")
open(f)
while(length(line <- readLines(f,n=1)) > 0) {
  pileup_string <- strsplit(line, "\t")[[1]][5]
  
  pileup_string <- gsub("n", "N", pileup_string)
  pileup_string <- gsub("a", "A", pileup_string)
  pileup_string <- gsub("t", "T", pileup_string)
  pileup_string <- gsub("c", "C", pileup_string)
  pileup_string <- gsub("g", "G", pileup_string)
  
  pileup_string <- gsub("#", "*", pileup_string)
  pileup_string <- gsub(",", ".", pileup_string)
  
  ref.seq <- unlist(strsplit(line, "\t")[[1]][3])
  
  ins <- insert.text(pileup_string)
  del <- del.text(pileup_string)
  del.ref.out <- del.ref(pileup_string, ref.seq)
  noInDel <- noIndel.text(pileup_string)
  
  this.dat <- table(strsplit(noInDel,"")[[1]])
  starts.out <- as.character(this.dat[names(this.dat)=="^"])
  if(length(starts.out)==0){starts.out=0}
  noInDelASCII <- gsub("(?:\\^).{1}", "", noInDel)
  
  ins.table <- table(ins)
  ins.out <- paste0(mapply(x=names(ins.table), y=ins.table, function(x,y) paste0(x,":",y)), collapse="|")
  
  del.table <- table(del)
  del.out <- paste0(mapply(x=names(del.table), y=del.table, function(x,y) paste0(x,":",y)), collapse="|")
  
  del.ref.table <- table(del.ref.out)
  del.ref.table.out <- paste0(mapply(x=names(del.ref.table), y=del.ref.table, function(x,y) paste0(x,":",y)), collapse="|")
  
  this.dat <- table(strsplit(noInDelASCII,"")[[1]])
  
  A.out=as.numeric(this.dat[names(this.dat)=="A"])
  if(ref.seq=="A"){A.out=this.dat[names(this.dat)=="."]}
  if(length(A.out)==0){A.out=0}
  
  T.out=as.numeric(this.dat[names(this.dat)=="T"])
  if(ref.seq=="T"){T.out=this.dat[names(this.dat)=="."]}
  if(length(T.out)==0){T.out=0}
  
  C.out=as.numeric(this.dat[names(this.dat)=="C"])
  if(ref.seq=="C"){C.out=this.dat[names(this.dat)=="."]}
  if(length(C.out)==0){C.out=0}
  
  G.out=as.numeric(this.dat[names(this.dat)=="G"])
  if(ref.seq=="G"){G.out=this.dat[names(this.dat)=="."]}
  if(length(G.out)==0){G.out=0}
  
  N.out=as.numeric(this.dat[names(this.dat)=="N"])
  if(ref.seq=="N"){N.out=this.dat[names(this.dat)=="."]}
  if(length(N.out)==0){N.out=0}
  
  stops.out <- as.character(this.dat[names(this.dat)=="$"])
  if(length(stops.out)==0){stops.out=0}
  
  out.table <- paste(unlist(strsplit(line, "\t")[[1]][1]),
                     unlist(strsplit(line, "\t")[[1]][2]),
                     unlist(strsplit(line, "\t")[[1]][3]),
                     unlist(strsplit(line, "\t")[[1]][4]),
                     A.out,
                     T.out,
                     C.out,
                     G.out,
                     N.out,
                     starts.out,
                     stops.out,
                     ins.out,
                     del.out,
                     del.ref.table.out, sep="\t")
  cat(out.table, sep = "\n")
  
}
