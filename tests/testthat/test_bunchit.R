# testing the main bunchit function

# testing that it works with kink, n_boot = 0
kink1 <- bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 100,
                 bins_l = 20, bins_r = 20, t0 = 0, t1 = 1, n_boot = 0)
test_that("bunchit returns correct objects with n_boot = 0", {
    expect_is(kink1, "list")
    expect_is(kink1$plot, "ggplot")
    expect_is(kink1$data, "data.frame")
    expect_is(kink1$model_fit, "matrix")
    expect_length(kink1, 22)
    expect_length(kink1$B, 1)
    expect_length(kink1$B_vector, 1)
    expect_length(kink1$B_sd, 1)
    expect_length(kink1$b, 1)
    expect_length(kink1$b_vector, 1)
    expect_length(kink1$b_sd, 1)
    expect_length(kink1$e, 1)
    expect_length(kink1$e_vector, 1)
    expect_length(kink1$e_sd, 1)
    expect_length(kink1$alpha, 1)
    expect_length(kink1$alpha_vector, 1)
    expect_length(kink1$alpha_sd, 1)
    expect_length(kink1$zD, 1)
    expect_length(kink1$zD_bin, 1)
    expect_length(kink1$zU_bin, 1)
    expect_length(kink1$marginal_buncher, 1)
    expect_length(kink1$marginal_buncher_vector, 1)
    expect_length(kink1$marginal_buncher_sd, 1)
})

# testing that it works with kink, n_boot = 100,
# returns vector for b_vector, B_vector, e_vector, marginal_buncher,
# and 100 NAs for alpha_vector
kink2 <- bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 100,
                 bins_l = 20, bins_r = 20, t0 = 0, t1 = 1, n_boot = 100)
test_that("bunchit returns correct objects with positive n_boot", {
    expect_length(kink2$B_vector, 100)
    expect_length(kink2$B_sd, 1)
    expect_length(kink2$b, 1)
    expect_length(kink2$b_vector, 100)
    expect_length(kink2$b_sd, 1)
    expect_length(kink2$e, 1)
    expect_length(kink2$e_vector, 100)
    expect_length(kink2$e_sd, 1)
    expect_length(kink2$alpha, 1)
    expect_length(kink2$alpha_vector, 100)
    expect_equal(sum(is.na(kink2$alpha_vector)), 100) # 100 NAs since no alpha
    expect_length(kink2$alpha_sd, 1)
    expect_length(kink2$marginal_buncher, 1)
    expect_length(kink2$marginal_buncher_vector, 100)
    expect_length(kink2$marginal_buncher_sd, 1)
})

# test that non-parametric elasticity is correct (kink)
kink3 <- bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 100,
                 bins_l = 20, bins_r = 20, t0 = 0, t1 = .2)
elasticity <- ((kink3$marginal_buncher - 10000)/10000)/.2
test_that("non-parametric elasticity in kink setting is correct", {
    expect_equal(kink3$e, round(elasticity,3))
})


# test that parametric elasticity is correct (kink)
kink4 <- bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 100,
                 bins_l = 20, bins_r = 20, t0 = 0, t1 = .2, e_parametric = T)
Dz <- kink4$b * 100
Dz_over_zstar <- Dz/10000
dt <- .2
t0 = 0
elasticity <- -log(1+Dz_over_zstar)/log(1-(dt/(1-t0)))
test_that("non-parametric elasticity in kink setting is correct", {
    expect_equal(kink4$e, round(elasticity,3))
})



# testing that it works with notch, n_boot = 0
notch1 <- bunchit(z_vector = bunching_data$notch_vector, zstar = 10000, binwidth = 50,
                  bins_l = 40, bins_r = 40, t0 = 0.18, t1 = .25, notch = TRUE, correct = FALSE, n_boot = 0)
test_that("bunchit returns correct objects with n_boot = 0", {
    expect_is(notch1, "list")
    expect_is(notch1$plot, "ggplot")
    expect_is(notch1$data, "data.frame")
    expect_is(notch1$model_fit, "matrix")
    expect_length(notch1, 22)
    expect_length(notch1$B, 1)
    expect_length(notch1$B_vector, 1)
    expect_length(notch1$B_sd, 1)
    expect_length(notch1$b, 1)
    expect_length(notch1$b_vector, 1)
    expect_length(notch1$b_sd, 1)
    expect_length(notch1$e, 1)
    expect_length(notch1$e_vector, 1)
    expect_length(notch1$e_sd, 1)
    expect_length(notch1$alpha, 1)
    expect_length(notch1$alpha_vector, 1)
    expect_length(notch1$alpha_sd, 1)
    expect_length(notch1$zD, 1)
    expect_length(notch1$zD_bin, 1)
    expect_length(notch1$zU_bin, 1)
    expect_length(notch1$marginal_buncher, 1)
    expect_length(notch1$marginal_buncher_vector, 1)
    expect_length(notch1$marginal_buncher_sd, 1)
})

# testing that it works with notch, n_boot > 0, (returns correct objects)
# and also gives correct warning that zD is above marginal buncher
# trigger this by setting correct = TRUE
test_that("returns all objects and correct warning that zD > marginal buncher (notch case)", {
    expect_warning(notch2 <- bunchit(z_vector = bunching_data$notch_vector, zstar = 10000, binwidth = 50,
                                     bins_l = 40, bins_r = 40, t0 = 0.18, t1 = .25, correct = TRUE, notch = TRUE))
    expect_length(notch2$B_vector, 100)
    expect_length(notch2$B_sd, 1)
    expect_length(notch2$b, 1)
    expect_length(notch2$b_vector, 100)
    expect_length(notch2$b_sd, 1)
    expect_length(notch2$e, 1)
    expect_length(notch2$e_vector, 100)
    expect_length(notch2$e_sd, 1)
    expect_length(notch2$alpha, 1)
    expect_length(notch2$alpha_vector, 100)
    expect_length(notch2$alpha_sd, 1)
    expect_length(notch2$marginal_buncher, 1)
    expect_length(notch2$marginal_buncher_vector, 100)
    expect_length(notch2$marginal_buncher_sd, 1)
})


# test that bunching with correction is smaller than bunching without
test_that("bunching with correction is smaller than without", {
    b1 <- bunchit(z_vector = bunching_data$notch_vector, zstar = 10000, binwidth = 50,
                  bins_l = 40, bins_r = 40, t0 = 0.18, t1 = .25, notch = TRUE,
                  correct = FALSE, n_boot = 0)$b
    b2 <- suppressWarnings(bunchit(z_vector = bunching_data$notch_vector, zstar = 10000, binwidth = 50,
                                   bins_l = 40, bins_r = 40, t0 = 0.18, t1 = .25, notch = TRUE,
                                   correct = TRUE, n_boot = 0)$b)
    expect_gt(b1,b2)
})

# test that bunching with correct_above_zu correction is smaller than without
test_that("bunching with correct_above_zu is smaller than without", {
    b1 <- suppressWarnings(bunchit(z_vector = bunching_data$notch_vector, zstar = 10000, binwidth = 50,
                                   bins_l = 40, bins_r = 40, t0 = 0.18, t1 = .25, notch = TRUE,
                                   correct = TRUE, correct_above_zu = FALSE, n_boot = 0)$b)
    b2 <- suppressWarnings(bunchit(z_vector = bunching_data$notch_vector, zstar = 10000, binwidth = 50,
                                   bins_l = 40, bins_r = 40, t0 = 0.18, t1 = .25, notch = TRUE,
                                   correct = TRUE, correct_above_zu = TRUE, n_boot = 0)$b)
    expect_gt(b1,b2)
})

# test that force_notch works
test_that("force_notch works", {
    zu_bin <- suppressWarnings(bunchit(z_vector = bunching_data$notch_vector, zstar = 10000, binwidth = 50,
                                       bins_l = 40, bins_r = 40, t0 = 0.18, t1 = .25, notch = TRUE,
                                       correct = FALSE, n_boot = 0, bins_excl_r = 2, force_notch = TRUE)$zU_bin)
    expect_equal(zu_bin, 2)
})


# test non-parametric elasticity for notches
notch3 <- bunchit(z_vector = bunching_data$notch_vector, zstar = 10000, binwidth = 50,
                  bins_l = 40, bins_r = 40, t0 = 0.18, t1 = 0.25, notch = TRUE,
                  correct = FALSE, n_boot = 0)
Dz <- notch3$b * 50
Dz_over_zstar <- Dz/10000
dt = 0.25 - 0.18
t0 = .18
elasticity <- (1/(2+Dz_over_zstar))*(Dz_over_zstar**2)/(dt/(1-t0))
test_that("non-parametric elasticity in notch setting is correct", {
    expect_equal(notch3$e, round(elasticity,3))
})

# test that parametric elasticity for notches fails correctly (hits lower bound)
test_that("parametric elasticity correctly hits lower bound when no solution", {
    notch4 <- suppressWarnings(bunchit(z_vector = bunching_data$notch_vector, zstar = 10000, binwidth = 50,
                                       bins_l = 40, bins_r = 40, t0 = 0.9, t1 = 0.95, notch = TRUE,
                                       correct = FALSE, n_boot = 0, e_parametric = TRUE))
    expect_equal(notch4$e, 0)
})


# test if invalid extra_fe correctly caught
test_that("invalid extra_fe correctly caught", {
    expect_warning(bunchit(z_vector = bunching_data$notch_vector, zstar = 10000, binwidth = 50,
                      bins_l = 40, bins_r = 40, t0 = 0.18, t1 = 0.25, notch = TRUE,
                      correct = FALSE, n_boot = 0, extra_fe = 11075),
                   "Warning: extra_fe do not correspond to any bin value. \n Consider using the plot_hist() function and check the returned data's bin values", fixed = TRUE)
})



# cleanup
rm(kink1, kink2, kink3, kink4, elasticity, Dz, Dz_over_zstar, dt, t0, notch1, notch3)
