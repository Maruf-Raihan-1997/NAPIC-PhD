# DH% (OPA assay) paired t-test

pre  <- c(0.207, 0.210, 0.170)
post <- c(0.396, 0.436, 0.457)

# Paired t-test
t.test(post, pre, paired = TRUE)

# Optional: compute differences
diff <- post - pre
mean(diff)
sd(diff)

# TPC (Folin assay) paired t-test

pre  <- c(22.12, 20.99, 23.15)
post <- c(80.04, 65.48, 64.15)

# Paired t-test
t.test(post, pre, paired = TRUE)

# Optional: compute differences
diff <- post - pre
mean(diff)
sd(diff)