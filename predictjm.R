is.factor.df <- function(x) {
  out <- as.vector(NULL)
  for(i in 1:ncol(x)) {
    out <- c(out, is.factor(x[,i]))
  }
  return(out)
}

predict.jm <- function(m) {
  m_dat <- attr(m, 'frame')
  m_main <- attr(terms(m),'term.labels')[!grepl(':',attr(terms(m),'term.labels'))]
  m_main_f <- m_main[is.factor.df(m_dat[,m_main])]
  m_main_c <- m_main[!is.factor.df(m_dat[,m_main])]
}

# predict.jm <- function(x, newdata=data.frame(NULL), method='0') {
#   x_dat <- attr(x,'frame')
#   x_main <- attr(terms(x),'term.labels')[!grepl(':',attr(terms(x),'term.labels'))]
#   x_dat_main <- subset(x_dat, select=x_main)
#   x_dat_factors <- as.data.frame(x_dat_main[,is.factor.df(x_dat_main)])
#   names(x_dat_factors) <- names(x_dat_main)[is.factor.df(x_dat_main)]
#   x_dat_covs <- x_dat_main[,!is.factor.df(x_dat_main)]
#   names(x_dat_covs) <- names(x_dat_main)[!is.factor.df(x_dat_main)]
#   
#   if(method=='0') {
#     
#   }
#   if(method!='0') stop('Method not implemented')
# }