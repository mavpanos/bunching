context("test_plot_hist_inputs")

# load the data
data("bunching_data")

# test that all required parameters have non-missing inputs

test_that("main parameter inputs are not missing", {
    # no input for z_vector
    expect_error(plot_hist(zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20))

    # no input for zstar
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, binwidth = 50,
                           bins_l = 20, bins_r = 20))

    # no input for binwidth
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000,
                           bins_l = 20, bins_r = 20))

    # no input for bins_l
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000,
                           binwidth = 50, bins_r = 20))

    # no input for bins_r
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000,
                           binwidth = 50, bins_r = 20))

})


# test that all parameters are correctly specified

test_that("main parameter inputs are correctly specified", {
    # z_vector is numeric
    expect_error(plot_hist(z_vector = c(1, NA), zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20))

    expect_error(plot_hist(z_vector = as.character(bunching_data$kink_vector), zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20))

    # binv is not min/max/median
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binv = "med", binwidth = 50,
                           bins_l = 20, bins_r = 20))

    # zstar is outside data
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 20000, binwidth = 50,
                           bins_l = 20, bins_r = 20),
                 "zstar is outside of z_vector's range of values")


    # binwidth is not positive
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = -10,
                           bins_l = 20, bins_r = 20),
                 "Binwidth must be a positive number")

    # bins_l not a a positive integer
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = -20, bins_r = 20),
                 "bins_l must be a positive integer")

    # bins_r not a a positive integer
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = -20),
                 "bins_r must be a positive integer")


    # p_title_size not positive
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_title_size = -5),
                 "p_title_size must be a positive numeric value")

    # p_title_size not numeric
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_title_size = "five"),
                 "p_title_size must be a positive numeric value")

    # p_axis_title_size not positive
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_axis_title_size = -5),
                 "p_axis_title_size must be a positive numeric value")

    # p_axis_title_size not numeric
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_axis_title_size = "five"),
                 "p_axis_title_size must be a positive numeric value")

    # p_axis_val_size not positive
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_axis_val_size = -5),
                 "p_axis_val_size must be a positive numeric value")

    # p_axis_val_size not numeric
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_axis_val_size = "five"),
                 "p_axis_val_size must be a positive numeric value")

    # p_miny not numeric
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_miny = "five"),
                 "plot's minimum y-value p_miny must be numeric")

    # p_maxy not numeric
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_maxy = "five"),
                 "plot's maximum y-value p_miny must be numeric")

    # p_maxy < p_miny
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_miny = 100, p_maxy = 50),
                 "p_maxy cannot be smaller than p_miny")

    # p_ybreaks contains non-numeric elements (if of length 1)
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_ybreaks = "100"),
                 "p_ybreaks must only contain numeric values")

    # p_ybreaks contains character (if length > 1)
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_ybreaks = c("100",200)),
                 "p_ybreaks must only contain numeric values")

    # p_ybreaks contains NA (if length > 1)
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_ybreaks = c(200, NA)),
                 "p_ybreaks must only contain numeric values")

    # p_freq_size not numeric
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_freq_size = "5"),
                 "p_freq_size must be numeric")

    # p_freq_msize not a positive number
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_freq_msize = -1),
                 "p_freq_msize must be a positive numeric value")

    # p_freq_msize not numeric
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_freq_msize = "1"),
                 "p_freq_msize must be a positive numeric value")

    # p_zstar_size not a positive number
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_zstar_size = -1),
                 "p_zstar_size must be a positive numeric value")

    # p_zstar not logical
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_zstar = 1),
                 "p_zstar can only be TRUE or FALSE", fixed = TRUE)

    # p_zstar_size not numeric
    expect_error(plot_hist(z_vector = bunching_data$kink_vector, zstar = 10000, binwidth = 50,
                           bins_l = 20, bins_r = 20, p_zstar_size = "1"),
                 "p_zstar_size must be a positive numeric value")
})

