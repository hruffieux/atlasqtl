library(ggplot2)

# Relationship between diff_lb and e:
adjusted_lognormal_cdf <- function(x, mu, sigma, m) {
  return(m + (1 - m) * pnorm((log(x) - mu) / sigma))
}


data.frame(
   diff_lb = seq(0, 50, length.out = 5000),
   e = adjusted_lognormal_cdf(seq(0, 50, length.out = 5000), mu=0.1, sigma=0.5, m=0)
) %>% ggplot(aes(x = diff_lb, y = e))+
  geom_point()+
  theme_bw()+
  scale_y_continuous(limits = c(0, 1))



# Relationship between iterations and e:
shifted_logic <- function(x, k, x_0, m) {
  m + (1 - m) * (1 - 1 / (1 + exp(-k * (x - x_0))))
}

data.frame(
  iteration = 1:500,
  e = shifted_logic(1:500, k=0.1, x_0=50, m = 0.1)
) %>% ggplot(aes(x = iteration, y = e))+
  geom_point()+
  theme_bw()+
  scale_y_continuous(limits = c(0, 1))

