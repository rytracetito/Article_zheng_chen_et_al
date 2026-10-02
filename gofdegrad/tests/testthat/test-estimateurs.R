test_that("estimer_wiener_homogene recupere les parametres theoriques", {
    set.seed(1)
    T <- 1:50
    th0 <- list(mu_W = 0.3, sigma = 0.08)
    X <- sim_proc(T, th0, "Wiener", mp = 20000)
    est <- estimer_wiener_homogene(X, T)
    expect_equal(est$mu_W, th0$mu_W, tolerance = 0.01)
    expect_equal(est$sigma, th0$sigma, tolerance = 0.01)
  })

test_that("estimer_gamma_homogene recupere les parametres theoriques", {
    set.seed(1)
    T <- 1:50
    th0 <- list(alpha = 1.5, beta = 0.7)
    X <- sim_proc(T, th0, "Gamma", mp = 20000)
    est <- estimer_gamma_homogene(X, T)
    expect_equal(est$alpha, th0$alpha, tolerance = 0.02)
    expect_equal(est$beta, th0$beta, tolerance = 0.02)
  })

test_that("estimer_ig_homogene recupere les parametres theoriques", {
    set.seed(1)
    T <- 1:50
    th0 <- list(mu_IG = 0.8, eta = 1.2)
    X <- sim_proc(T, th0, "IG", mp = 20000)
    est <- estimer_ig_homogene(X, T)
    expect_equal(est$mu_IG, th0$mu_IG, tolerance = 0.01)
    expect_equal(est$eta, th0$eta, tolerance = 0.02)
  })

test_that("estimer_wiener_nonhomog recupere les parametres theoriques", {
    set.seed(1)
    T <- 1:50
    th0 <- list(mu_W = 0.15, sigma = 0.4)
    X <- sim_proc(T, th0, "Wiener", mp = 20000, Lam = function(t) t^1.7)
                    est <- estimer_wiener_nonhomog(X, T)
                    expect_equal(est$mu_W, th0$mu_W, tolerance = 0.01)
                    expect_equal(est$sigma, th0$sigma, tolerance = 0.01)
                    expect_equal(est$r, 1.7, tolerance = 0.01)
                  })

  test_that("estimer_gamma_nonhomog recupere les parametres theoriques", {
      set.seed(1)
      T <- 1:50
      th0 <- list(alpha = 0.8, beta = 0.5)
      X <- sim_proc(T, th0, "Gamma", mp = 20000, Lam = function(t) t^1.4)
                      est <- estimer_gamma_nonhomog(X, T)
                      expect_equal(est$alpha, th0$alpha, tolerance = 0.02)
                      expect_equal(est$beta, th0$beta, tolerance = 0.02)
                      expect_equal(est$r, 1.4, tolerance = 0.02)
                    })

    test_that("estimer_parametres route vers le bon estimateur", {
        set.seed(1)
        T <- 1:50
        X <- sim_proc(T, list(mu_W = 0.2, sigma = 0.05), "Wiener", mp = 6)
        expect_identical(
              estimer_parametres(X, T, "Wiener"),
              estimer_wiener_homogene(X, T)
            )
      })

    test_that("estimer_parametres refuse un modele inconnu", {
        T <- 1:10
        X <- matrix(1:10, 1, 10)
        expect_error(estimer_parametres(X, T, "inexistant"), "unknown model")
      })

    test_that("les estimateurs Gamma/IG detectent des increments negatifs", {
        T <- 1:10
        X <- matrix(c(1, 0.5, 0.3, rep(1, 7)), nrow = 1)  # increment negatif en position 2
        expect_error(estimer_gamma_homogene(X, T), "negative increments")
        expect_error(estimer_ig_homogene(X, T), "negative increments")
      })

    test_that("increments_positifs tolere un increment nul avec avertissement", {
        T <- 1:5
        X <- matrix(c(1, 1, 2, 3, 4), nrow = 1)  # increment nul en position 2 (1 -> 1)
        expect_warning(dX <- gofdegrad:::increments_positifs(X, "test", "Gamma"), "zero increment")
        expect_true(all(dX > 0))
      })
    
