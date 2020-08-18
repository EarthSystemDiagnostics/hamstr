
BlockAverage <- function(ts1, scale){

  ts2 <- stats::filter(ts1, rep(1/scale, scale))
  ts2 <- ts(ts2[seq(ceiling(scale/2), length(ts2), by = scale)])
  return(ts2)
}

GetARCoef <- function(ts1){
  f1 <- arima0(ts1, order = c(1,0,0))
  return(f1$coef[1])
}


ar1 <- 0.8
scale <- 1:20
ts1 <-  arima.sim(list(ar = ar1), 1e04)
#plot(ts1)
ar1.coefs <- sapply(scale, function(x){
  ts2 <- BlockAverage(ts1, x)
  GetARCoef(ts2)
})

plot(scale, ar1.coefs, ylim = c(-0.2,1), type = "l")
abline(h = c(ar1, 0))
lines(scale, (ar1^scale), col = "Blue")



library(tidyverse)

set.seed(20200818)
n <- 30
dat <- tibble(
  x = 1:n,
  y = as.numeric(arima.sim(list(ar = 0.9), n)),
  p = ifelse(x %in% seq(2, length(x), by = 3), TRUE, FALSE)
)

dat2 <- dat %>%
  do({
    tibble(
      x = seq(2, length(.$y), 3),
      y = BlockAverage(.$y, 3)
    )
  })

dat %>%
  ggplot(aes(x=x, y=y)) +
  geom_line(aes(colour = "All points")) +
  geom_point(aes(colour = "All points"))+
  scale_colour_manual("", values = c("All points" = "Darkgrey"),
                        drop = FALSE) +
  theme_bw()


dat %>%
  ggplot(aes(x=x, y=y)) +
  #geom_line() +
  geom_point(aes(colour = "All points")) +
  geom_line(aes(colour = "All points")) +
  geom_point(data = filter(dat, p == TRUE), aes(colour = "Lag 3"))+
  geom_line(data = filter(dat, p == TRUE), aes(colour = "Lag 3"))+
  geom_point(data = dat2, aes(colour = "Block Averaged"))+
  geom_line(data = dat2, aes(colour = "Block Averaged")) +
  theme_bw() +
  scale_colour_manual("", values = c("All points" = "Grey", "Block Averaged" = "Orange", "Lag 3" = "Blue")) +
  geom_vline(xintercept = seq(0.5, n, by = 3), colour = "Lightgrey") +
  theme(panel.grid = element_blank())




BA.AR1 <- crossing(rep = 1:10, AR1 = seq(0.5, 0.9, by = 0.2),
                scl = 1:20) %>%
  group_by(rep, AR1, scl) %>%
  do({
    ts1 <- arima.sim(list(ar = .$AR1), 1e04)
    tibble(
      AR1.naive = .$AR1^.$scl,
      AR1.BlockAverage = GetARCoef(BlockAverage(ts1, scale = .$scl[1]))
    )
  })

fig.BA.AR1 <- BA.AR1 %>%
  ungroup() %>%
  gather(type, AR1.coef, -scl, -AR1, -rep) %>%
  group_by(type, AR1, scl) %>%
  summarise_if(is.numeric, mean) %>%
  #filter(AR1 == 0.1) %>%
  ggplot(aes(x = scl, y = AR1.coef, colour = type, group = paste(type, rep))) +
  geom_line() +
  facet_wrap(~AR1, ncol = 1, labeller = label_both) +
  labs(x = "Block width") +
  expand_limits(x  = 0)

ggsave("doc/fig-BA-AR1.png", fig.BA.AR1, height = 4, width = 4.5)
