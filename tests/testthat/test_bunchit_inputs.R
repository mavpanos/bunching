context("test_bunchit_inputs")

# load the data
data("bunching_data")

# test that all required parameters have non-missing inputs

test_that("main parameter inputs are not missing", {
    # no input for z_vector
    expect_error(bunchit(zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2))

    # no input for zstar
    expect_error(bunchit(z_vector = bunching_data$kink_vector, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2))

    # no input for binwidth
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2))

    # no input for bins_l
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000,
                         binwidth = 50, bins_r = 20, t0 = 0, t1 = 0.2))

    # no input for bins_r
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000,
                         binwidth = 50, bins_r = 20, t0 = 0, t1 = 0.2))

    # no input for t0
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000,
                         binwidth = 50, bins_l = 20, bins_r = 20, t1 = 0.2))

    # no input for t1
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000,
                         binwidth = 50, bins_l = 20, bins_r = 20, t0 = 0.2))

})


# test that all parameters are correctly specified

test_that("main parameter inputs are correctly specified", {
    # z_vector is numeric
    expect_error(bunchit(z_vector = c(1, NA), zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2))

    expect_error(bunchit(z_vector = as.character(bunching_data$kink_vector), zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2))

    # binv is not min/max/median
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binv = "med", binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2))

    # zstar is outside data
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 20000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2),
                 "zstar is outside of z_vector's range of values")

    # zstar  is 0
    expect_error(bunchit(z_vector = bunching_data$kink_vector - 10000, zstar = 0, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2),
                 "zstar cannot be zero. If this your true bunching point, must re-centre it away from zero")

    # binwidth is not positive
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = -10,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2),
                 "Binwidth must be a positive number")

    # bins_l not a a positive integer
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = -20, bins_r = 20, t0 = 0, t1 = 0.2),
                 "bins_l must be a positive integer")

    # bins_r not a a positive integer
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = -20, t0 = 0, t1 = 0.2),
                 "bins_r must be a positive integer")

    # poly must be non-negative
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, poly = - 5),
                 "poly must be a non-negative integer")

    # poly must be integer
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, poly = 5.5),
                 "poly must be a non-negative integer")

    # excluded bins below zstar cannot be negative
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, bins_excl_l = -2),
                 "Number of bins in bunching region below zstar must be a non-negative integer")

    # excluded bins below zstar must be an integer
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, bins_excl_l = 2.5),
                 "Number of bins in bunching region below zstar must be a non-negative integer")

    # excluded bins above zstar cannot be negative
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, bins_excl_r = -2),
                 "Number of bins in bunching region above zstar must be a non-negative integer")

    # excluded bins above zstar must be an integer
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, bins_excl_r = 2.5),
                 "Number of bins in bunching region above zstar must be a non-negative integer")

    # bunching region below zstar too wide
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, bins_excl_l = 15),
        "Bunching region below zstar too wide. Increase bins_l relative to bins_excl_l")

    # bunching region above zstar too wide
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, bins_excl_r = 15),
                 "Bunching region above zstar too wide. Increase bins_r relative to bins_excl_r")

    # rn includes zero
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, rn = c(0,50)),
                 "rn cannot include zero as a round number")

    # rn includes non-round numbers
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, rn = c(20.5,50)),
                 "Round number(s) must be integer(s)", fixed = TRUE) # fixed = T required for strings with brackets

    # rn includes more than two numbers
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, rn = c(25,50,100)),
                 "rn cannot include more than two unique levels for round number bunching")

    # rn includes two identical numbers
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, rn = c(100,100)),
                 "the two round numbers in rn cannot be identical")

    # rn includes numbers ouside of data range
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, rn = c(100,20000)),
                 "rn includes round numbers outside of z_vector's range of values")

    # rn includes numbers that are too large for z_vector's range
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, rn = 6000),
                 "rn includes round numbers that are too large for z_vector's range of values. \n Use smaller values or increase bins_l and bins_r")

    # number of bootstrap iterations is not a non-negative integer
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, n_boot = -1),
                 "bootstrap sample size n_boot must be a non-negative integer")

    # number of bootstrap iterations is not a non-negative integer
    expect_warning(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, n_boot = 50),
                 "Are you sure your choice for the bootstrap sample size n_boot is large enough?", fixed = TRUE)

    # correct (integration constraint) not logical
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, correct = 1),
                 "correct (for integration constraint) can only be TRUE or FALSE", fixed = TRUE)

    # correct_above_zu not logical
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, correct_above_zu = 1),
                  "correct_above_zu (for integration constraint in notch setting) can only be TRUE or FALSE", fixed = TRUE)

    # correct_iter_max not positive
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = 0.2, correct_iter_max = 0),
                 "Maximum number of iterations (for integration constraint) correct_iter_max must be a positive integer", fixed = TRUE)

    # t0 not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = "zero", t1 = 0.2),
                 "t0 must be a numeric value")

    # t1 not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = 0, t1 = "one"),
                 "t1 must be a numeric value")

    # t0 unusually low
    expect_warning(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = -0.1, t1 = .5),
                   "Are you sure this is the correct value for t0? \n Note that t0 = 0 means a 0% and t0 = 1 a 100% tax rate!", fixed = TRUE)

    # t0 unusually high
    expect_warning(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, t0 = 1.1, t1 = 1.5),
                   "Are you sure this is the correct value for t0? \n Note that t0 = 0 means a 0% and t0 = 1 a 100% tax rate!", fixed = TRUE)

    # t1 unusually low
    expect_warning(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, t0 = -1, t1 = -0.5),
                   "Are you sure this is the correct value for t0? \n Note that t0 = 0 means a 0% and t0 = 1 a 100% tax rate!", fixed = TRUE)

    # t1 unusually high
    expect_warning(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, t0 = 1.1, t1 = 1.5),
                   "Are you sure this is the correct value for t1? \n Note that t1 = 0 means a 0% and t1 = 1 a 100% tax rate!", fixed = TRUE)

    # t0 = t1
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, t0 = .1, t1 = .1),
                   "Cannot calculate elasticity (t0 cannot equal t1", fixed = TRUE)

    # t0 > t1
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, t0 = .1, t1 = .05),
                 "t1 must be larger than t0")

    # notch choice not logical
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, notch = 1),
                 "notch can either be TRUE or FALSE (i.e. kink)", fixed = TRUE)

    # force_notch choice not logical
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, force_notch = 1),
                 "force_notch can either be TRUE or FALSE")

    # p_title_size not positive
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_title_size = -5),
                 "p_title_size must be a positive numeric value")

    # p_title_size not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_title_size = "five"),
                 "p_title_size must be a positive numeric value")

    # p_axis_title_size not positive
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_axis_title_size = -5),
                 "p_axis_title_size must be a positive numeric value")

    # p_axis_title_size not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_axis_title_size = "five"),
                 "p_axis_title_size must be a positive numeric value")

    # p_axis_val_size not positive
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_axis_val_size = -5),
                 "p_axis_val_size must be a positive numeric value")

    # p_axis_val_size not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_axis_val_size = "five"),
                 "p_axis_val_size must be a positive numeric value")

    # p_miny not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_miny = "five"),
                 "plot's minimum y-value p_miny must be numeric")

    # p_maxy not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_maxy = "five"),
                 "plot's maximum y-value p_miny must be numeric")

    # p_maxy < p_miny
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_miny = 100, p_maxy = 50),
                 "p_maxy cannot be smaller than p_miny")

    # p_ybreaks contains non-numeric elements (if of length 1)
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_ybreaks = "100"),
                 "p_ybreaks must only contain numeric values")

    # p_ybreaks contains character (if length > 1)
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_ybreaks = c("100",200)),
                 "p_ybreaks must only contain numeric values")

    # p_ybreaks contains NA (if length > 1)
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_ybreaks = c(200, NA)),
                 "p_ybreaks must only contain numeric values")

    # p_freq_size not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_freq_size = "5"),
                 "p_freq_size must be numeric")

    # p_freq_msize not a positive number
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_freq_msize = -1),
                 "p_freq_msize must be a positive numeric value")

    # p_freq_msize not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_freq_msize = "1"),
                 "p_freq_msize must be a positive numeric value")

    # p_cf_size not a positive number
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_cf_size = -1),
                 "p_cf_size must be a positive numeric value")

    # p_cf_size not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_cf_size = "1"),
                 "p_cf_size must be a positive numeric value")

    # p_zstar_size not a positive number
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_zstar_size = -1),
                 "p_zstar_size must be a positive numeric value")

    # p_zstar_size not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_zstar_size = "1"),
                 "p_zstar_size must be a positive numeric value")

    # p_b not logical
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_b = 1),
                 "p_b (whether to show bunching estimate on plot) must be TRUE or FALSE", fixed = TRUE)

    # p_b_e_xpos not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_b_e_xpos = TRUE),
                 "p_b_e_xpos must be numeric and lie within z_vector's range")

    # p_b_e_xpos not within z_vector's range
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_b_e_xpos = 20000),
                 "p_b_e_xpos must be numeric and lie within z_vector's range")

    # p_b_e_ypos not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_b_e_ypos = TRUE),
                 "p_b_e_ypos must be numeric")

    # p_b_e_size not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_b_e_size = TRUE),
                 "p_b_e_size must be a positive numeric value")

    # p_b_e_size not positive
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, p_b_e_size = -10),
                 "p_b_e_size must be a positive numeric value")

    # seed not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, seed = TRUE),
                 "seed must be numeric")

    # e_parametric not logical
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, e_parametric = 1),
                 "e_parametric can only be TRUE or FALSE")

    # e_parametric_lb not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, e_parametric_lb = TRUE),
                 "e_parametric_lb must be numeric")

    # e_parametric_ub not numeric
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2, e_parametric_ub = TRUE),
                 "e_parametric_ub must be numeric")

    # e_parametric_lb > e_parametric_ub
    expect_error(bunchit(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                         bins_l = 20, bins_r = 20, t0 = .1, t1 = .2,
                         e_parametric_lb = 1, e_parametric_ub = .5),
                 "e_parametric_ub must be larger than e_parametric_lb")
    })

