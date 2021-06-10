## s
test_that("Model matrix for a single population curve", {
  tm <- get_transform_bs(4, 2)                  # transform matrix
  dmi <- get_design_bs(1:4, 2, 1)               # design matrix info
  sobj <- s(1:4, knots = 2, degree = 1)
  expect_mapequal(sobj, list(degree = 1,
                             type = 'bs',
                             is_sub = FALSE,
                             knots = dmi$knots,
                             trans_mat = tm[, -1],
                             model_mat = structure(dmi$design %*% tm[, -1],
                                                   penalty = cbind(0, diag(2)),
                                                   spl_dim = 3,
                                                   is_sub = FALSE,
                                                   level = '1')))
})


test_that("Model matrix for multiple population curves", {
  f <- as.factor(c(1, 1, 1, 1, 2, 2, 2, 2))
  tm <- get_transform_bs(8, 4)                  # transform matrix
  dmi <- get_design_bs(c(1:4, 1:4), 4, 3)               # design matrix info
  sobj <- s(c(1:4, 1:4), by = f, knots = 4, degree = 3)
  mm <- Matrix::.bdiag(split.data.frame(dmi$design %*% tm[, -1], f))
  mm <- methods::as(mm, 'dgCMatrix')
  expect_mapequal(sobj, list(degree = 3,
                             type = 'bs',
                             is_sub = FALSE,
                             knots = dmi$knots,
                             trans_mat = tm[, -1],
                             model_mat = structure(mm,
                                                   penalty = cbind(0, 0, 0, diag(4)),
                                                   spl_dim = 7,
                                                   is_sub = FALSE,
                                                   level = c('1', '2'))))
})

test_that("Model matrix for subject-specific curves", {
  f <- as.factor(c(1, 1, 1, 1, 2, 2, 2, 2))
  tm <- get_transform_bs(5, 3)     # transform matrix
  dmi <- get_design_bs(c(1:4, 2:5), 2, 2) # design matrix info
  sobj <- s(c(1:4, 2:5), by = f, knots = 2, degree = 2, is_sub = TRUE)
  expect_mapequal(sobj, list(degree = 2,
                             type = 'bs',
                             is_sub = TRUE,
                             knots = dmi$knots,
                             trans_mat = tm,
                             model_mat = structure(split.data.frame(dmi$design %*% tm, f),
                                                   block_dim = 3,
                                                   spl_dim = 5,
                                                   is_sub = TRUE,
                                                   level = c('1', '2'))))
})




## parse_spline




## parse_effect



## parse_response




