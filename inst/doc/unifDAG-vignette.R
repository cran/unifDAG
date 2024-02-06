## ----FMversion, echo = FALSE--------------------------------------------------
UDversion <- packageDescription('unifDAG')$Version

## ----echo=FALSE,warning=FALSE,message=FALSE,cache=FALSE-------
options(digits=7)
options(width=64)
require('knitr')

library(unifDAG)

opts_chunk$set(
   keep.source=FALSE,
   out.width='4.5in',
   fig.width=10,
   fig.height=6,
   fig.path = 'figure/',
   cache=FALSE,
   tidy=TRUE,
   fig.align='center',
   size='small')

## ----unifDAGExpl, echo = TRUE---------------------------------
myWgtFun <- function(m,lB,uB) { runif(m,lB,uB) }
set.seed(123)
dag1 <- unifDAG(n = 10, weighted = TRUE, wFUN = list(myWgtFun, 0, 2))
dag2 <- unifDAG.approx(n = 150, n.exact = 40, weighted = TRUE,
                       wFUN = list(myWgtFun, 0, 2))

## ----simpleExpl, echo = TRUE----------------------------------
cnt <- c(0,0,0) ## count occurances of G1, G2 and G3
set.seed(123)
for (i in 1:100) {
  g <- unifDAG(n = 2, weighted = FALSE)
  m <- as(g, "matrix") ## adjacency matrix
  if ( (m[2,1]==0) & (m[1,2]==0) ) {
    cnt[3] <- cnt[3] + 1 ## G3
  } else if (m[2,1]==0) {
    cnt[1] <- cnt[1] + 1 ## G1
  } else {
    cnt[2] <- cnt[2] + 1 ## G2
  }
}
cnt

