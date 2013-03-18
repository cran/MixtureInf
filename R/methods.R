print.emtest <- function(x, digits=max(3, getOption('digits') - 3), ...) {
  family_label <- switch(x$family,
                         normal='Normal mixture with known variance 1',
                         normal_unequalvar='Normal mixture without equal variance assumption',
                         poisson='Poisson mixture',
                         paste(x$family, 'mixture'))
    
  cat('\nEM-test for H_0: order=', x$m0, 'under', family_label, '\n\n')

  cat('MLE of mixing distribution under null hypothesis:\n')
  cat('Mixing proportions:', format(x$alpha, digits), '\n')
  cat('Mixing parameters:', format(x$theta, digits), '\n\n')

  cat('EM-test statistics:', format(x$emstat, digits), '\n')
  cat('p-values:', format(x$pvalue, digits), '\n')

  invisible(x)
}
