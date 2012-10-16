getLoglike.mat <- function(y, x1.num, n) {

    ret <- matrix(data=FALSE, nrow=n, ncol=6)

    #if (gmodel3) {

      col <- 3*y + 1 + x1.num
      #col <- 3*D + 1 + snp

      #> cbind(y,x1.num)[1:8,]
      #     y x1.num
      #[1,] 0      1
      #[2,] 0      1
      #[3,] 0      1
      #[4,] 0      1
      #[5,] 0      0
      #[6,] 0      2
      #[7,] 0      2         #         1      2      3       4      5     6
      #[8,] 1      1      --> (D,G)= (0,0), (0,1), (0,2), (1,0), (1,1), (1,2)

      #> col[1:10]
      # [1] 2 2 2 2 1 3 3 5 2 3     ok
    #}

    #else {
    #  col <- 3*D + 1 + fsnp
    #}

    for (i in 1:n) ret[i, col[i]] <- TRUE
    ret

} # END: getLoglike.mat


#> tt=getLoglike.mat(y, x1.num, n)
#      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]
#[1,] FALSE  TRUE FALSE FALSE FALSE FALSE   ---> this means first row individual has D=1 and G=1
#[2,] FALSE  TRUE FALSE FALSE FALSE FALSE
#[3,] FALSE  TRUE FALSE FALSE FALSE FALSE
#[4,] FALSE  TRUE FALSE FALSE FALSE FALSE
#[5,]  TRUE FALSE FALSE FALSE FALSE FALSE
#> cbind(y,x1.num)[1:5,]
#     y x1.num
#[1,] 0      1
#[2,] 0      1
#[3,] 0      1
#[4,] 0      1
#[5,] 0      0          -->  column 1 is TRUE in the fifth row above
#