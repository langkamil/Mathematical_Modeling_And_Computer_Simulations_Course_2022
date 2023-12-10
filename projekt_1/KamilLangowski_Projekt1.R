# Projekt 1
# Aleksandra Burczyk
# Monika Gmaj
# Kinga Jankowska
# Kamil Langowski


# Cz. I -------------------------------------------------------------------




# Zadanie 1.  -------------------------------------------------------------


# Ciagi liczb pseudolosowych generowane metoda dystrybuanty odwrotnej dla 
# zadanej gestosci (rozklad jednostajny U(0,1))


pacman::p_load(moments, pracma)
# Rozklad jednostajny
set.seed(1)
sim_unif <- runif(10^3,0,1)



# 1.1 wyznaczyc w sposob analityczny dystrybuante i dystrybuante odwrotna

################# ARCUS SINUS #####################

# rozklad Arcus Sinus F^(x) : U(0,1) -> R
Arcsin <- function(x){ 2/pi*asin(sqrt(x))}


# rozklad Arcus Sinus F^-1(x) : U(0,1) -> R
inv_Arcsin <- function(x){ 1/2 - 1/2*cos(x*pi) }

Dane_Arcsin1 <- inv_Arcsin(runif(10^2,0,1))
Dane_Arcsin2 <- inv_Arcsin(runif(10^3,0,1))
Dane_Arcsin3 <- inv_Arcsin(runif(10^5,0,1))


# 1.2 przedstawienie dystrybuanty na wykresie
dystrybuanty <- function(x1,x2,d1,d2, ylimits) {
  plot(x1,d1(x1),type = "l",ylim = ylimits,col="green",ylab = "", lwd=2,
       main = "Funkcja gestosci i gestosci odwrotnej")
  lines(x2,d2(x2), lwd=2)
  legend("topleft",inset=.02, legend=c("F(x)", "F^-1(x)"),
         col=c("green", "black"),lty=3, cex=0.75,
         box.lty=1)
}

par(mfrow=c(1,1))
n = 1000
dystrybuanty(seq(0,1,length=n), seq(0,1,length=n), Arcsin, inv_Arcsin, c(0,1))

# 1.3 Miary statystyczne 
Miary <- function(Dane, t_params){
  Miary_Statystyczne <- c("Mediana", "Srednia", "Dominanta", "Kurtoza",
                          "Wariancja", "Odchylenie standardowe", "Rozstep")
  Wyniki_Z_Generatora <- c(median(Dane), mean(Dane),density(Dane)$x[which.max(density(Dane)$y)],
                           kurtosis(Dane), var(Dane), std(Dane), max(Dane)-min(Dane))
  Wyniki_Teoretyczne <- t_params
  tabela <- data.frame(Miary_Statystyczne,Wyniki_Z_Generatora,Wyniki_Teoretyczne)
  print(summary(Dane))
  print(tabela)
}


Miary(Dane_Arcsin1, c(1/2, 1/2, '{0,1}', 3/2, 1/8, sqrt(1/8), 1))
Miary(Dane_Arcsin2, c(1/2, 1/2, '{0,1}', 3/2, 1/8, sqrt(1/8), 1))
Miary(Dane_Arcsin3, c(1/2, 1/2, '{0,1}', 3/2, 1/8, sqrt(1/8), 1))



# 1.4 Histogram i gestosc
hist_PDF <- function(dane){
  par(mfrow=c(1,2))
  hist(dane, prob=TRUE, breaks=40,xlab = "x", main='Histogram',ylab = "Czestotliwosc")
  plot(density(dane), type = "l", pch = 19, 
       col = "red", xlab = "x", ylab = "Gestosc", main='Gestość')
}

hist_PDF(Dane_Arcsin1)
hist_PDF(Dane_Arcsin2)
hist_PDF(Dane_Arcsin3)

# 1.5 Dystrybuanta i Dystrubuanta empiryczna
dyst_emp_teor <- function(d, d_inv, n, xlims, umin=0, umax=1, tmin=0, tmax=1) {
  plot(ecdf(d_inv(runif(n,umin,umax))),do.points=FALSE,
       main=paste("Dystrybuanta empiryczna i zgodna z rozkladem"),
       ylab = "", xlim = xlims)
  
  x <- seq(tmin,tmax,0.01)
  lines(x,d(x),type = "l",col="red")
  legend("bottomright",inset=.02, legend=c("F(x)-Dyst. emp", "F(x)-Dyst rozkl"),
         col=c("black", "red"),lty=3, cex=0.75,
         box.lty=1)
}

par(mfrow=c(1,1))
dyst_emp_teor(Arcsin,inv_Arcsin, n=10^2, xlims = c(0,1))
dyst_emp_teor(Arcsin,inv_Arcsin, n=10^3, xlims = c(0,1))
dyst_emp_teor(Arcsin,inv_Arcsin, n=10^5, xlims = c(0,1))

################# CAUCHY #####################

# rozklad Cauchyego F^(x) : U(0,1) -> R

Cauchy <- function(x,x0 = 0 ,alfa = 1) 1/2 + (1/pi)*atan((x-x0)/alfa)

# rozklad Cauchyego F^-1(x) : U(0,1) -> R

inv_Cauchy <- function(x,x0 = 0,alfa = 1/2){ x0 - alfa*cot(pi*x)}


Dane_inv_Cauchy <- inv_Cauchy(sim_unif)
Dane_inv_Cauchy <- Dane_inv_Cauchy[Dane_inv_Cauchy > -10 & Dane_inv_Cauchy < 10] # usunac mocno odstajace obs

Dane_Cauchy1 <- inv_Cauchy(runif(10^2,0,1))
Dane_Cauchy2 <- inv_Cauchy(runif(10^3,0,1))
Dane_Cauchy3 <- inv_Cauchy(runif(10^5,0,1))



# 1.2 przedstawienie dystrybuanty na wykresie
par(mfrow = c(1,1))
dystrybuanty(seq(-5,5,length=n),seq(0,1,length=n), Cauchy, inv_Cauchy, c(-3,3))

# 1.3 Miary statystyczne 


Miary(Dane_Cauchy1, c(0,0,0,'nieoznaczony','nieoznaczony','nieoznaczony','infty'))
Miary(Dane_Cauchy2, c(0,0,0,'nieoznaczony','nieoznaczony','nieoznaczony','infty'))
Miary(Dane_Cauchy3, c(0,0,0,'nieoznaczony','nieoznaczony','nieoznaczony','infty'))


# 1.4 Histogram i gestosc

hist_PDF(Dane_Cauchy1)
hist_PDF(Dane_Cauchy2)
hist_PDF(Dane_Cauchy3)

# 1.5 Dystrybuanta i Dystrubuanta empiryczna

par(mfrow=c(1,1))
dyst_emp_teor(Cauchy,inv_Cauchy, n=1000, xlims= c(-10,10), umin=-10, umax=10, tmin=-10, tmax=10)


################# WARTOSCI EKSTREMALNYCH #####################

# rozklad wartosci ekstremalnych F^(x) : U(0,1) -> R
extr.val <- function(x,alfa=1,beta=1/2) 1-exp(-exp(x-alfa)/beta)

# rozklad wartosci ekstremalnych F^-1(x) : U(0,1) -> R
inv_extr.val <- function(x,alfa=1,beta=1/2) {(log(-alfa*log((1-x)))/beta)}

hist(inv_extr.val(sim_unif,alfa=1,beta=1/2), prob=TRUE, breaks=40)
lines(density(inv_extr.val(sim_unif,alfa=1,beta=1/2)))


Dane_extr.val1 <- inv_extr.val(runif(10^2,0,1))
Dane_extr.val2 <- inv_extr.val(runif(10^3,0,1))
Dane_extr.val3 <- inv_extr.val(runif(10^5,0,1))

# 1.2 przedstawienie dystrybuanty na wykresie
dystrybuanty(seq(-5,5,length=n),seq(0,1,length=n), extr.val, inv_extr.val, c(-2,3))

# 1.3 Miary statystyczne 


Miary(Dane_extr.val1, c(0,1-1/2*-digamma(1),1,0,1/2^2*pi^2/6,0,0))
Miary(Dane_extr.val2, c(0,1-1/2*-digamma(1),1,0,1/2^2*pi^2/6,0,0))
Miary(Dane_extr.val3, c(0,1-1/2*-digamma(1),1,0,1/2^2*pi^2/6,0,0))

# 1.4 Histogram i gestosc
hist_PDF(Dane_extr.val1)
hist_PDF(Dane_extr.val2)
hist_PDF(Dane_extr.val3)

# 1.5 Dystrybuanta i Dystrybuanta empiryczna

par(mfrow=c(1,1))
dyst_emp_teor(extr.val,inv_extr.val, n=10000, xlims= c(-5,5), umin=-5, umax=5, tmin=-5, tmax=5)



################# Rayleigh #####################
# rozklad Reyleigha F^(x) : U(0,1) -> R
Rayleigh <- function(x,alfa = 0,beta=1)  1-exp(-((x-alfa)/beta)^2) 


# rozklad Reyleigha F^-1(x) : U(0,1) -> R
inv_Rayleigh <- function(x,alfa = 0,beta=1)  alfa+beta*(sqrt(log(1/(1-x))))

hist(inv_Rayleigh(sim_unif), prob=TRUE, breaks=40)
lines(density(inv_Rayleigh(sim_unif)))

Dane_Rayleigh1 <- inv_Rayleigh(runif(10^2,0,1))
Dane_Rayleigh2 <- inv_Rayleigh(runif(10^3,0,1))
Dane_Rayleigh3 <- inv_Rayleigh(runif(10^5,0,1))

# 1.2 przedstawienie dystrybuanty na wykresie
dystrybuanty(seq(0,1,length=n),seq(0,1,length=n), Rayleigh, inv_Rayleigh, c(-2,3))

# 1.3 Miary statystyczne 



Miary(Dane_Rayleigh1, c(0+1*sqrt(log(2)),0+1*sqrt(pi)/2,0+1/sqrt(2),(32-3*pi^2)/(4-pi)^2,1^2*(1-pi/4),sqrt(1^2*(1-pi/4)),'infty'))
Miary(Dane_Rayleigh2, c(0+1*sqrt(log(2)),0+1*sqrt(pi)/2,0+1/sqrt(2),(32-3*pi^2)/(4-pi)^2,1^2*(1-pi/4),sqrt(1^2*(1-pi/4)),'infty'))
Miary(Dane_Rayleigh3, c(0+1*sqrt(log(2)),0+1*sqrt(pi)/2,0+1/sqrt(2),(32-3*pi^2)/(4-pi)^2,1^2*(1-pi/4),sqrt(1^2*(1-pi/4)),'infty'))

# 1.4 Histogram i gestosc

hist_PDF(Dane_Rayleigh1)
hist_PDF(Dane_Rayleigh2)
hist_PDF(Dane_Rayleigh3)

# 1.5 Dystrybuanta i Dystrubuanta empiryczna

par(mfrow=c(1,1))
dyst_emp_teor(Rayleigh,inv_Rayleigh, n=1000, xlims= c(0,3), umin=0, umax=1 , tmax=3)


################# Weibull #####################

# rozklad Weibulla F^(x) : U(0,1) -> R
Weibull <- function(x,alfa=1,beta=2,c=1)  1-exp(-((x-alfa)/beta)^c) # x > alfa

# rozklad Weibulla F^-1(x) : U(0,1) -> R
inv_Weibull <- function(x,alfa=1,beta=2,c=1) alfa + beta*((-log(1-x))^(1/c))


Dane_Weibull1 <- inv_Weibull(runif(10^2,0,1))
Dane_Weibull2 <- inv_Weibull(runif(10^3,0,1))
Dane_Weibull3 <- inv_Weibull(runif(10^5,0,1))

# 1.2 przedstawienie dystrybuanty na wykresie
dystrybuanty(seq(0,1,length=n),seq(0,1,length=n), Weibull, inv_Weibull, c(-3,10))

# 1.3 Miary statystyczne 


Miary(Dane_Weibull1, c(1+2*(log(2)^(1/1)),1+2*gamma((1+1)/1),1,(2^2/1)*(2*gamma(2)-(gamma(1))^2),(2^2)*(gamma(3)-(gamma(2))^2),sqrt((2^2)*(gamma(3)-(gamma(2))^2)),'infty'))
Miary(Dane_Weibull2, c(1+2*(log(2)^(1/1)),1+2*gamma((1+1)/1),1,(2^2/1)*(2*gamma(2)-(gamma(1))^2),(2^2)*(gamma(3)-(gamma(2))^2),sqrt((2^2)*(gamma(3)-(gamma(2))^2)),'infty'))
Miary(Dane_Weibull3, c(1+2*(log(2)^(1/1)),1+2*gamma((1+1)/1),1,(2^2/1)*(2*gamma(2)-(gamma(1))^2),(2^2)*(gamma(3)-(gamma(2))^2),sqrt((2^2)*(gamma(3)-(gamma(2))^2)),'infty'))

# 1.4 Histogram i gestosc


hist_PDF(Dane_Weibull1)
hist_PDF(Dane_Weibull2)
hist_PDF(Dane_Weibull3)

# 1.5 Dystrybuanta i Dystrubuanta empiryczna 

par(mfrow=c(1,1))
dyst_emp_teor(Weibull,inv_Weibull, n=100, xlims= c(0,10), umin=0, umax=1, tmin=1, tmax=10)



################# LOGISTYCZNY #####################

# rozklad Logistyczny F^(x) : U(0,1) -> R
Logist <- function(x,alfa = 0 ,beta = 1) 1/(1+exp(-(x-alfa)/beta))

# rozklad Logistyczny F^-1(x) : U(0,1) -> R

inv_Logist <- function(x,alfa = 0,beta = 1) alfa + beta*log(x/(1-x))

hist(inv_Logist(sim_unif), prob=TRUE, breaks=40)
lines(density(inv_Logist(sim_unif)))

Dane_Logist1 <- inv_Logist(runif(10^2,0,1))
Dane_Logist2 <- inv_Logist(runif(10^3,0,1))
Dane_Logist3 <- inv_Logist(runif(10^5,0,1))

# 1.2 przedstawienie dystrybuanty na wykresie
dystrybuanty(seq(-5,5,length=n),seq(0,1,length=n), Logist, inv_Logist, c(-3,3))

# 1.3 Miary statystyczne 



Miary(Dane_Logist1, c(0,0,0,pi^2*1/3,pi^2/3*1^2,sqrt(pi^2/3*1^2),'infty'))
Miary(Dane_Logist2, c(0,0,0,pi^2*1/3,pi^2/3*1^2,sqrt(pi^2/3*1^2),'infty'))
Miary(Dane_Logist3, c(0,0,0,pi^2*1/3,pi^2/3*1^2,sqrt(pi^2/3*1^2),'infty'))

# 1.4 Histogram i gestosc


hist_PDF(Dane_Logist1)
hist_PDF(Dane_Logist2)
hist_PDF(Dane_Logist3)

# 1.5 Dystrybuanta i Dystrubuanta empiryczna 

par(mfrow=c(1,1))
dyst_emp_teor(Logist,inv_Logist, n=100, xlims= c(-5,7), umin=0, umax=1, tmin=-5, tmax=7)


################# PARETO #####################

# rozklad Pareto F^(x) : U(0,1) -> R
Pareto <- function(x,c = 3)  1-x^(-c)

# rozklad Pareto F^-1(x) : U(0,1) -> R
inv_Pareto <- function(x,c = 3) (1-x)^(-1/c)

hist(inv_Pareto(sim_unif), prob=TRUE, breaks=40)
lines(density(inv_Pareto(sim_unif)))

Dane_Pareto1 <- inv_Pareto(runif(10^2,0,1))
Dane_Pareto2 <- inv_Pareto(runif(10^3,0,1))
Dane_Pareto3 <- inv_Pareto(runif(10^5,0,1))

# 1.2 przedstawienie dystrybuanty na wykresie
par(mfrow=c(1,1))
dystrybuanty(seq(0,1,length=n),seq(0,1,length=n),Pareto, inv_Pareto, c(-3,3))

# 1.3 Miary statystyczne 


Miary(Dane_Pareto1, c(2^(1/3),3/2,1,0,3-(3/2)^2,sqrt(3-(3/2)^2),'infty'))
Miary(Dane_Pareto2, c(2^(1/3),3/2,1,0,3-(3/2)^2,sqrt(3-(3/2)^2),'infty'))
Miary(Dane_Pareto3, c(2^(1/3),3/2,1,0,3-(3/2)^2,sqrt(3-(3/2)^2),'infty'))

# 1.4 Histogram i gestosc

hist_PDF(Dane_Pareto1)
hist_PDF(Dane_Pareto2)
hist_PDF(Dane_Pareto3)

# 1.5 Dystrybuanta i Dystrubuanta empiryczna 

par(mfrow=c(1,1))
dyst_emp_teor(Pareto,inv_Pareto, n=1000, xlims= c(0,50), umin=0, umax=1, tmin=-5, tmax=50)



################# TROJKATNY #####################

# rozklad Pareto F^(x) : U(0,1) -> R

Triangular_base <- function(x,a = 0,c = 1/2, b= 1) {
  if (x < a) value <- 0
  else if (x <= c) value <-(x-a)^2/((b-a)*(c-a))
  else if (x < b) value <- 1 - (b-x)^2/((b-a)*(b-c))
  else if (b <= x) value <- 1
  
  return(value)
}


Triangular <- function(x, ...) {
  sapply(x, Triangular_base, ...)
}


# rozklad Trojkatny F^-1(x) : U(0,1) -> R

inv_Triangular_base <- function(x,a = 0,c = 1/2, b= 1) {
  if (x < (c - a) / (b - a)) value <- a + sqrt(max(0, (c - a) * (b - a) * x))
  else if (x >= (c - a) / (b - a)) value <- b - sqrt(max(0, (b - a) * (b - c) * (1 - x)))
  return(value)
} 
inv_Triangular <- function(x, ...) {
  sapply(x, inv_Triangular_base, ...)
}

hist(inv_Triangular(sim_unif), prob=TRUE, breaks=40)
lines(density(inv_Triangular(sim_unif)))

Dane_Triangular1 <- inv_Triangular(runif(10^2,0,1))
Dane_Triangular2 <- inv_Triangular(runif(10^3,0,1))
Dane_Triangular3 <- inv_Triangular(runif(10^5,0,1))
# 1.2 przedstawienie dystrybuanty na wykresie
dystrybuanty(seq(0,1,length=n),seq(0,1,length=n), Triangular, inv_Triangular, c(0,1))

# 1.3 Miary statystyczne 


Miary(Dane_Triangular1, c(min(Dane_Triangular3)+sqrt((max(Dane_Triangular3)-min(Dane_Triangular3))*(1/2 - min(Dane_Triangular3))/2),(min(Dane_Triangular3)+max(Dane_Triangular3)+1/2)/3,1/2,-3/5,(3*(max(Dane_Triangular3)-min(Dane_Triangular3))^2+max(Dane_Triangular3)+min(Dane_Triangular3)-1)^2/72,sqrt((3*(max(Dane_Triangular3)-min(Dane_Triangular3))^2+max(Dane_Triangular3)+min(Dane_Triangular3)-1)^2/72),max(Dane_Triangular3)-min(Dane_Triangular3))) 
Miary(Dane_Triangular2, c(min(Dane_Triangular3)+sqrt((max(Dane_Triangular3)-min(Dane_Triangular3))*(1/2 - min(Dane_Triangular3))/2),(min(Dane_Triangular3)+max(Dane_Triangular3)+1/2)/3,1/2,-3/5,(3*(max(Dane_Triangular3)-min(Dane_Triangular3))^2+max(Dane_Triangular3)+min(Dane_Triangular3)-1)^2/72,sqrt((3*(max(Dane_Triangular3)-min(Dane_Triangular3))^2+max(Dane_Triangular3)+min(Dane_Triangular3)-1)^2/72),max(Dane_Triangular3)-min(Dane_Triangular3)))
Miary(Dane_Triangular3, c(min(Dane_Triangular3)+sqrt((max(Dane_Triangular3)-min(Dane_Triangular3))*(1/2 - min(Dane_Triangular3))/2),(min(Dane_Triangular3)+max(Dane_Triangular3)+1/2)/3,1/2,-3/5,(3*(max(Dane_Triangular3)-min(Dane_Triangular3))^2+max(Dane_Triangular3)+min(Dane_Triangular3)-1)^2/72,sqrt((3*(max(Dane_Triangular3)-min(Dane_Triangular3))^2+max(Dane_Triangular3)+min(Dane_Triangular3)-1)^2/72),max(Dane_Triangular3)-min(Dane_Triangular3)))

# 1.4 Histogram i gestosc

hist_PDF(Dane_Triangular1)
hist_PDF(Dane_Triangular2)
hist_PDF(Dane_Triangular3)

# 1.5 Dystrybuanta i Dystrubuanta empiryczna 

par(mfrow=c(1,1))
dyst_emp_teor(Triangular,inv_Triangular, n=1000, xlims= c(0,1), umin=0, umax=1, tmin=0, tmax=1)



# Zadanie 2 ---------------------------------------------------------------

# Dla wygenerowanych próbek z Zadania 1 ciągów liczb pseudolosowych o wybranym rozkładzie
# 2.1 Znormalizować próbki

normalizacja <- function(x) (x - min(x))/(max(x) - min(x))

{
  Dane_Arcsin1_normalizowane <- normalizacja(Dane_Arcsin1)
  Dane_Arcsin2_normalizowane <- normalizacja(Dane_Arcsin2)
  Dane_Arcsin3_normalizowane <- normalizacja(Dane_Arcsin3)
}

{
  Dane_Rayleigh1_normalizowane <- normalizacja(Dane_Rayleigh1)
  Dane_Rayleigh2_normalizowane <- normalizacja(Dane_Rayleigh2)
  Dane_Rayleigh3_normalizowane <- normalizacja(Dane_Rayleigh3)
}




# 2.2
# Dla każdej znormalizowanej próbki wygenerować ciąg liczb pseudolosowych U wykorzystując (2),
# gdzie fN oznacza N-tą iterację mapy chaotycznej (1) dla p = 0.45.

# definiujemy funkcje mapy chaotycznej z polecenia
ukosna_mapa_namiotu <- function(zk, p=0.45) {
  if(zk >= 0 & zk < p) {
    return(zk/p)
  }
  if(zk >= p & zk <= 1) {
    return((1-zk)/(1-p))
  }
  else{
    print('Nieprawidłowa wartość')
  }
}


# definiujemy N-te zlozenie mapy chaotycznej z polecelnia
ukosna_mapa_namiotu_N_iteracji <- function(N, x0, p = 0.45){
  x = c()
  for (j in 1:length(x0)) {
    x[j] = ukosna_mapa_namiotu(x0[j], p)
    
    if(N > 1) {
      for(i in 2:N){
        x[j] = ukosna_mapa_namiotu(x[j], p)
      }
    }
  }
  return(x)
}



{
  Dane_Arcsin1_prawdziwie_losowe_1 <- ukosna_mapa_namiotu_N_iteracji(1, Dane_Arcsin1_normalizowane)
  Dane_Arcsin2_prawdziwie_losowe_1 <- ukosna_mapa_namiotu_N_iteracji(1, Dane_Arcsin2_normalizowane)
  Dane_Arcsin3_prawdziwie_losowe_1 <- ukosna_mapa_namiotu_N_iteracji(1, Dane_Arcsin3_normalizowane)
}

{
  Dane_Arcsin1_prawdziwie_losowe_5 <- ukosna_mapa_namiotu_N_iteracji(5, Dane_Arcsin1_normalizowane)
  Dane_Arcsin2_prawdziwie_losowe_5 <- ukosna_mapa_namiotu_N_iteracji(5, Dane_Arcsin2_normalizowane)
  Dane_Arcsin3_prawdziwie_losowe_5 <- ukosna_mapa_namiotu_N_iteracji(5, Dane_Arcsin3_normalizowane)
}

{
  Dane_Arcsin1_prawdziwie_losowe_12 <- ukosna_mapa_namiotu_N_iteracji(12, Dane_Arcsin1_normalizowane)
  Dane_Arcsin2_prawdziwie_losowe_12 <- ukosna_mapa_namiotu_N_iteracji(12, Dane_Arcsin2_normalizowane)
  Dane_Arcsin3_prawdziwie_losowe_12 <- ukosna_mapa_namiotu_N_iteracji(12, Dane_Arcsin3_normalizowane)
}


{
  Dane_Rayleigh1_prawdziwie_losowe_1 <- ukosna_mapa_namiotu_N_iteracji(1, Dane_Rayleigh1_normalizowane)
  Dane_Rayleigh2_prawdziwie_losowe_1 <- ukosna_mapa_namiotu_N_iteracji(1, Dane_Rayleigh2_normalizowane)
  Dane_Rayleigh3_prawdziwie_losowe_1 <- ukosna_mapa_namiotu_N_iteracji(1, Dane_Rayleigh3_normalizowane)
}

{
  Dane_Rayleigh1_prawdziwie_losowe_5 <- ukosna_mapa_namiotu_N_iteracji(5, Dane_Rayleigh1_normalizowane)
  Dane_Rayleigh2_prawdziwie_losowe_5 <- ukosna_mapa_namiotu_N_iteracji(5, Dane_Rayleigh2_normalizowane)
  Dane_Rayleigh3_prawdziwie_losowe_5 <- ukosna_mapa_namiotu_N_iteracji(5, Dane_Rayleigh3_normalizowane)
}


{
  Dane_Rayleigh1_prawdziwie_losowe_12 <- ukosna_mapa_namiotu_N_iteracji(12, Dane_Rayleigh1_normalizowane)
  Dane_Rayleigh2_prawdziwie_losowe_12 <- ukosna_mapa_namiotu_N_iteracji(12, Dane_Rayleigh2_normalizowane)
  Dane_Rayleigh3_prawdziwie_losowe_12 <- ukosna_mapa_namiotu_N_iteracji(12, Dane_Rayleigh3_normalizowane)
}

# 2.3 Wygenerować histogramy rozkładu U dla każdej z próbek, dla trzech wartości N = {1, 5, 12}

hist(Dane_Arcsin1_prawdziwie_losowe_1)
hist(Dane_Arcsin2_prawdziwie_losowe_1)
hist(Dane_Arcsin3_prawdziwie_losowe_1)

hist(Dane_Arcsin1_prawdziwie_losowe_5)
hist(Dane_Arcsin2_prawdziwie_losowe_5)
hist(Dane_Arcsin3_prawdziwie_losowe_5)

hist(Dane_Arcsin1_prawdziwie_losowe_12)
hist(Dane_Arcsin2_prawdziwie_losowe_12)
hist(Dane_Arcsin3_prawdziwie_losowe_12)




hist(Dane_Rayleigh1_prawdziwie_losowe_1)
hist(Dane_Rayleigh2_prawdziwie_losowe_1)
hist(Dane_Rayleigh3_prawdziwie_losowe_1)

hist(Dane_Rayleigh1_prawdziwie_losowe_5)
hist(Dane_Rayleigh2_prawdziwie_losowe_5)
hist(Dane_Rayleigh3_prawdziwie_losowe_5)

hist(Dane_Rayleigh1_prawdziwie_losowe_12)
hist(Dane_Rayleigh2_prawdziwie_losowe_12)
hist(Dane_Rayleigh3_prawdziwie_losowe_12)


# 2.4 W każdym przypadku wygenerować wykresy zależności par (uk,uk+1)
# (rozrzut punktów na płaszczyźnie (uk,uk+1)).

wykres_rozrzutu <- function(df) plot(x = df[-length(df)], y = df[-1],
                                     xlab='U_{k}',
                                     ylab = 'U_{k-1}',
                                     main = 'Wykres rozrzutu')

wykres_rozrzutu(Dane_Arcsin2_normalizowane)
wykres_rozrzutu(Dane_Arcsin2_prawdziwie_losowe_1)
wykres_rozrzutu(Dane_Arcsin2_prawdziwie_losowe_5)
wykres_rozrzutu(Dane_Arcsin2_prawdziwie_losowe_12)



wykres_rozrzutu(Dane_Rayleigh2_normalizowane)
wykres_rozrzutu(Dane_Rayleigh2_prawdziwie_losowe_1)
wykres_rozrzutu(Dane_Rayleigh2_prawdziwie_losowe_5)
wykres_rozrzutu(Dane_Rayleigh2_prawdziwie_losowe_12)





# Zadanie 3 ---------------------------------------------------------------
# Dla wygenerowanych próbek z Zadania 1 ciągów liczb pseudolosowych o wybranym rozkładzie
# wykonać podpunkty 2.1-2.4 dla dwóch dowolnie wybranych z [1] map chaotycznych.

#
# definiujemy funkcje mapy logistycznej
mapa_logistyczna <- function(x0, r=3.8) return(r*x0*(1 - x0))


# definiujemy N-te zlozenie mapy chaotycznej 
mapa_logistyczna_N_iteracji <- function(N, x0, r=3.8){
  x = c()
  for (j in 1:length(x0)) {
    
    x[j] = mapa_logistyczna(x0[j], r)
    if(N > 1){
      for(i in 2:N){
        x[j] = mapa_logistyczna(x[j], r)
      }
    }
  }
  return(x)
}


# definiujemy funkcje mapy kwadratowej
mapa_kwadratowa <- function(xk, a=4, b=0.5) return(b - a*xk^2)


# definiujemy N-te zlozenie mapy chaotycznej 
mapa_kwadratowa_N_iteracji <- function(N, x0, a=4, b=0.5){
  x = c()
  for (j in 1:length(x0)) {
    
    x[j] = mapa_kwadratowa(x0[j], a, b)
    if(N>1){
      for(i in 2:N){
        x[j] = mapa_kwadratowa(x[j], a, b)
      }
    }
  }
  return(x)
}

################# MAPA LOGISTYCZNA #####################


# 3.2 
# Dla każdej znormalizowanej próbki wygenerować ciąg liczb pseudolosowych U wykorzystując (2),
# gdzie fN oznacza N-tą iterację mapy chaotycznej (1) dla p = 0.45.

{
  Dane_Arcsin1_prawdziwie_losowe_1_logistyczna <- mapa_logistyczna_N_iteracji(1, Dane_Arcsin1_normalizowane)
  Dane_Arcsin2_prawdziwie_losowe_1_logistyczna <- mapa_logistyczna_N_iteracji(1, Dane_Arcsin2_normalizowane)
  Dane_Arcsin3_prawdziwie_losowe_1_logistyczna <- mapa_logistyczna_N_iteracji(1, Dane_Arcsin3_normalizowane)
}

{
  Dane_Arcsin1_prawdziwie_losowe_5_logistyczna <- mapa_logistyczna_N_iteracji(5, Dane_Arcsin1_normalizowane)
  Dane_Arcsin2_prawdziwie_losowe_5_logistyczna <- mapa_logistyczna_N_iteracji(5, Dane_Arcsin2_normalizowane)
  Dane_Arcsin3_prawdziwie_losowe_5_logistyczna <- mapa_logistyczna_N_iteracji(5, Dane_Arcsin3_normalizowane)
}

{
  Dane_Arcsin1_prawdziwie_losowe_12_logistyczna <- mapa_logistyczna_N_iteracji(12, Dane_Arcsin1_normalizowane)
  Dane_Arcsin2_prawdziwie_losowe_12_logistyczna <- mapa_logistyczna_N_iteracji(12, Dane_Arcsin2_normalizowane)
  Dane_Arcsin3_prawdziwie_losowe_12_logistyczna <- mapa_logistyczna_N_iteracji(12, Dane_Arcsin3_normalizowane)
}


{
  Dane_Rayleigh1_prawdziwie_losowe_1_logistyczna <- mapa_logistyczna_N_iteracji(1, Dane_Rayleigh1_normalizowane)
  Dane_Rayleigh2_prawdziwie_losowe_1_logistyczna <- mapa_logistyczna_N_iteracji(1, Dane_Rayleigh2_normalizowane)
  Dane_Rayleigh3_prawdziwie_losowe_1_logistyczna <- mapa_logistyczna_N_iteracji(1, Dane_Rayleigh3_normalizowane)
}

{
  Dane_Rayleigh1_prawdziwie_losowe_5_logistyczna <- mapa_logistyczna_N_iteracji(5, Dane_Rayleigh1_normalizowane)
  Dane_Rayleigh2_prawdziwie_losowe_5_logistyczna <- mapa_logistyczna_N_iteracji(5, Dane_Rayleigh2_normalizowane)
  Dane_Rayleigh3_prawdziwie_losowe_5_logistyczna <- mapa_logistyczna_N_iteracji(5, Dane_Rayleigh3_normalizowane)
}

{
  Dane_Rayleigh1_prawdziwie_losowe_12_logistyczna <- mapa_logistyczna_N_iteracji(12, Dane_Rayleigh1_normalizowane)
  Dane_Rayleigh2_prawdziwie_losowe_12_logistyczna <- mapa_logistyczna_N_iteracji(12, Dane_Rayleigh2_normalizowane)
  Dane_Rayleigh3_prawdziwie_losowe_12_logistyczna <- mapa_logistyczna_N_iteracji(12, Dane_Rayleigh3_normalizowane)
}


# 3.3 Wygenerować histogramy rozkładu U dla każdej z próbek, dla trzech wartości N = {1, 5, 12}

hist(Dane_Arcsin1_prawdziwie_losowe_1_logistyczna)
hist(Dane_Arcsin2_prawdziwie_losowe_1_logistyczna)
hist(Dane_Arcsin3_prawdziwie_losowe_1_logistyczna)

hist(Dane_Arcsin1_prawdziwie_losowe_5_logistyczna)
hist(Dane_Arcsin2_prawdziwie_losowe_5_logistyczna)
hist(Dane_Arcsin3_prawdziwie_losowe_5_logistyczna)

hist(Dane_Arcsin1_prawdziwie_losowe_12_logistyczna)
hist(Dane_Arcsin2_prawdziwie_losowe_12_logistyczna)
hist(Dane_Arcsin3_prawdziwie_losowe_12_logistyczna)




hist(Dane_Rayleigh1_prawdziwie_losowe_1_logistyczna)
hist(Dane_Rayleigh2_prawdziwie_losowe_1_logistyczna)
hist(Dane_Rayleigh3_prawdziwie_losowe_1_logistyczna)

hist(Dane_Rayleigh1_prawdziwie_losowe_5_logistyczna)
hist(Dane_Rayleigh2_prawdziwie_losowe_5_logistyczna)
hist(Dane_Rayleigh3_prawdziwie_losowe_5_logistyczna)

hist(Dane_Rayleigh1_prawdziwie_losowe_12_logistyczna)
hist(Dane_Rayleigh2_prawdziwie_losowe_12_logistyczna)
hist(Dane_Rayleigh3_prawdziwie_losowe_12_logistyczna)

# 3.4 W każdym przypadku wygenerować wykresy zależności par (uk,uk+1)
# (rozrzut punktów na płaszczyźnie (uk,uk+1)).


wykres_rozrzutu(Dane_Arcsin2_normalizowane)
wykres_rozrzutu(Dane_Arcsin2_prawdziwie_losowe_1_logistyczna)
wykres_rozrzutu(Dane_Arcsin2_prawdziwie_losowe_5_logistyczna)
wykres_rozrzutu(Dane_Arcsin2_prawdziwie_losowe_12_logistyczna)



wykres_rozrzutu(Dane_Rayleigh2_normalizowane)
wykres_rozrzutu(Dane_Rayleigh2_prawdziwie_losowe_1_logistyczna)
wykres_rozrzutu(Dane_Rayleigh2_prawdziwie_losowe_5_logistyczna)
wykres_rozrzutu(Dane_Rayleigh2_prawdziwie_losowe_12_logistyczna)



################# MAPA KWADRATOWA #####################

# 3.2 
# Dla każdej znormalizowanej próbki wygenerować ciąg liczb pseudolosowych U wykorzystując (2),
# gdzie fN oznacza N-tą iterację mapy chaotycznej (1) dla p = 0.45.

{
  Dane_Arcsin1_prawdziwie_losowe_1_kwadratowa <- mapa_kwadratowa_N_iteracji(1, Dane_Arcsin1_normalizowane)
  Dane_Arcsin2_prawdziwie_losowe_1_kwadratowa <- mapa_kwadratowa_N_iteracji(1, Dane_Arcsin2_normalizowane)
  Dane_Arcsin3_prawdziwie_losowe_1_kwadratowa <- mapa_kwadratowa_N_iteracji(1, Dane_Arcsin3_normalizowane)
}

{
  Dane_Arcsin1_prawdziwie_losowe_5_kwadratowa <- mapa_kwadratowa_N_iteracji(5, Dane_Arcsin1_normalizowane)
  Dane_Arcsin2_prawdziwie_losowe_5_kwadratowa <- mapa_kwadratowa_N_iteracji(5, Dane_Arcsin2_normalizowane)
  Dane_Arcsin3_prawdziwie_losowe_5_kwadratowa <- mapa_kwadratowa_N_iteracji(5, Dane_Arcsin3_normalizowane)
}

{
  Dane_Arcsin1_prawdziwie_losowe_12_kwadratowa <- mapa_kwadratowa_N_iteracji(12, Dane_Arcsin1_normalizowane)
  Dane_Arcsin2_prawdziwie_losowe_12_kwadratowa <- mapa_kwadratowa_N_iteracji(12, Dane_Arcsin2_normalizowane)
  Dane_Arcsin3_prawdziwie_losowe_12_kwadratowa <- mapa_kwadratowa_N_iteracji(12, Dane_Arcsin3_normalizowane)
}


{
  Dane_Rayleigh1_prawdziwie_losowe_1_kwadratowa <- mapa_kwadratowa_N_iteracji(1, Dane_Rayleigh1_normalizowane)
  Dane_Rayleigh2_prawdziwie_losowe_1_kwadratowa <- mapa_kwadratowa_N_iteracji(1, Dane_Rayleigh2_normalizowane)
  Dane_Rayleigh3_prawdziwie_losowe_1_kwadratowa <- mapa_kwadratowa_N_iteracji(1, Dane_Rayleigh3_normalizowane)
}

{
  Dane_Rayleigh1_prawdziwie_losowe_5_kwadratowa <- mapa_kwadratowa_N_iteracji(5, Dane_Rayleigh1_normalizowane)
  Dane_Rayleigh2_prawdziwie_losowe_5_kwadratowa <- mapa_kwadratowa_N_iteracji(5, Dane_Rayleigh2_normalizowane)
  Dane_Rayleigh3_prawdziwie_losowe_5_kwadratowa <- mapa_kwadratowa_N_iteracji(5, Dane_Rayleigh3_normalizowane)
}

{
  Dane_Rayleigh1_prawdziwie_losowe_12_kwadratowa <- mapa_kwadratowa_N_iteracji(12, Dane_Rayleigh1_normalizowane)
  Dane_Rayleigh2_prawdziwie_losowe_12_kwadratowa <- mapa_kwadratowa_N_iteracji(12, Dane_Rayleigh2_normalizowane)
  Dane_Rayleigh3_prawdziwie_losowe_12_kwadratowa <- mapa_kwadratowa_N_iteracji(12, Dane_Rayleigh3_normalizowane)
}


# 3.3 Wygenerować histogramy rozkładu U dla każdej z próbek, dla trzech wartości N = {1, 5, 12}

hist(Dane_Arcsin1_prawdziwie_losowe_1_kwadratowa)
hist(Dane_Arcsin2_prawdziwie_losowe_1_kwadratowa)
hist(Dane_Arcsin3_prawdziwie_losowe_1_kwadratowa)

hist(Dane_Arcsin1_prawdziwie_losowe_5_kwadratowa)
hist(Dane_Arcsin2_prawdziwie_losowe_5_kwadratowa)
hist(Dane_Arcsin3_prawdziwie_losowe_5_kwadratowa)

hist(Dane_Arcsin1_prawdziwie_losowe_12_kwadratowa)
hist(Dane_Arcsin1_prawdziwie_losowe_12_kwadratowa)
hist(Dane_Arcsin1_prawdziwie_losowe_12_kwadratowa)



hist(Dane_Rayleigh1_prawdziwie_losowe_1_kwadratowa)
hist(Dane_Rayleigh2_prawdziwie_losowe_1_kwadratowa)
hist(Dane_Rayleigh3_prawdziwie_losowe_1_kwadratowa)

hist(Dane_Rayleigh1_prawdziwie_losowe_5_kwadratowa)
hist(Dane_Rayleigh2_prawdziwie_losowe_5_kwadratowa)
hist(Dane_Rayleigh3_prawdziwie_losowe_5_kwadratowa)

hist(Dane_Rayleigh1_prawdziwie_losowe_12_kwadratowa)
hist(Dane_Rayleigh2_prawdziwie_losowe_12_kwadratowa)
hist(Dane_Rayleigh3_prawdziwie_losowe_12_kwadratowa)

# 3.4 W każdym przypadku wygenerować wykresy zależności par (uk,uk+1)
# (rozrzut punktów na płaszczyźnie (uk,uk+1)).


wykres_rozrzutu(Dane_Arcsin2_normalizowane)
wykres_rozrzutu(Dane_Arcsin2_prawdziwie_losowe_1_kwadratowa)
wykres_rozrzutu(Dane_Arcsin2_prawdziwie_losowe_5_kwadratowa)
wykres_rozrzutu(Dane_Arcsin2_prawdziwie_losowe_12_kwadratowa)



wykres_rozrzutu(Dane_Rayleigh2_normalizowane)
wykres_rozrzutu(Dane_Rayleigh2_prawdziwie_losowe_1_kwadratowa)
wykres_rozrzutu(Dane_Rayleigh2_prawdziwie_losowe_5_kwadratowa)
wykres_rozrzutu(Dane_Rayleigh2_prawdziwie_losowe_12_kwadratowa)

# Cz. II ------------------------------------------------------------------

library(OpenImageR)

# 1. Zapisujemy i ładujemy obrazek
img <- readImage("papugi_res.png")
imageShow(img)
# 2. Wyodrębniamy trzy macierze {R, G,B}

R <- img[,,1]
G <- img[,,2]
B <- img[,,3]
R_w <- c(R)
G_w <- c(G)
B_w <- c(B)

# 3. Dla każdego koloru oryginalnego obrazka możemy zobrazować intensywność 
# oraz przedstawić histogramy

hist(R_w, col = 'red', breaks = 150)
hist(G_w, col = 'green', breaks = 150)
hist(B_w, col = 'blue', breaks = 150)

plot(R_w*255, col = 'red', ylim = c(0,300), type = 'h')
plot(G_w*255, col = 'green', ylim = c(0,300), type = 'h')
plot(B_w*255, col = 'blue', ylim = c(0,300), type = 'h')

# 4. Wykorzystując mapę logistyczną dla wybranego r ∈ [3.6, 4] i x0 ∈ (0, 1) 
# konstruujemy wektor wartości o długości K = 3 ·N ·M. 
# Otrzymany wektor dzielimy na trzy podwektory {xR, xG, xB}. 
# Klucz kodowania składa się z dwóch informacji (x0, r).

# przyjmujemy x0 = 0.5 i r = 3.8

x0 <- 0.5
r <- 3.8
mapa_logistyczna <- function(x0, r=3.8) return(r*x0*(1 - x0))
mapa_logistyczna_wektor <- function(x, x0, r) {
  
  tmp <- numeric(length(x))
  
  tmp[1] <- x0
  
  for (i in 2:length(x)) {
    
    tmp[i] <- mapa_logistyczna(tmp[i-1], r)
    
  }
  
  return(tmp)
}

x_R <- mapa_logistyczna_wektor(R_w, x0 = 0.5, r =3.8)
x_G <- mapa_logistyczna_wektor(G_w, x0 = 0.5, r =3.8)
x_B <- mapa_logistyczna_wektor(B_w, x0 = 0.5, r =3.8)


# 5. Sortujemy (malejąco lub rosnąco) osobno wektory {xR, xG, xB} otrzymując {xsR, xsG, xsB}.
# Następnie tworzymy trzy wektory pozycji {PR, PG, PB}, które przechowują numery pozycji
# elementów z {xsR, xsG, xsB} w wektorach {xR, xG, xB}.




# sortowanie
x_R_s <- sort(x_R)
x_G_s <- sort(x_G)
x_B_s <- sort(x_B)




# wektory pozycji
P_R <- match(x_R_s, x_R)
P_G <- match(x_G_s, x_G)
P_B <- match(x_B_s, x_B)

# 6. Kodowanie
kodowanie <- function(x, P){
  
  tmp <- numeric(length(x))
  
  for (i in 1:length(x)){
    x1 <- x[i]
    x2 <- P[i]
    tmp[x2] <- x1
  }
  
  return(tmp)
}
# przygotowanie wektorów tasowanych

T_R <- kodowanie(R_w, P_R)
T_G <- kodowanie(G_w, P_G)
T_B <- kodowanie(B_w, P_B)

# 7. Łączenie wektorów

T <- array(data=c(T_R, T_G, T_B), dim = dim(img))

# wyświetlenie zakodowanego obrazka
imageShow(T)

# 8. Wykreślamy intensywność każdego koloru zakodowanego obrazka

plot(T_R*255, col = 'red', ylim = c(0,300), type = 'h')
plot(T_G*255, col = 'green', ylim = c(0,300), type = 'h')
plot(T_B*255, col = 'blue', ylim = c(0,300), type = 'h')


# 9. Dekodowanie

# klucz: x0 = 0.5 i r = 3.8


dekodowanie <- function(x, x0, r){
  x_mapowane <- mapa_logistyczna_wektor(x, x0, r)
  x_mapowane_sort <- sort(x_mapowane)
  x_mapowane_sort_pozycja <- match(x_mapowane_sort, x_mapowane)
  tmp <- numeric(length(x))
  for(i in (1:length(x))){
    pozycja <- x_mapowane_sort_pozycja[i]
    wartosc <- x[pozycja]
    tmp[i] <- wartosc
  }
  return(tmp)
}


D_R <- dekodowanie(T_R, 0.5, 3.8)
D_G <- dekodowanie(T_G, 0.5, 3.8)
D_B <- dekodowanie(T_B, 0.5, 3.8)

D <- array(c(D_R, D_G, D_B), dim=dim(img))


imageShow(D)


# 10  Wyznaczamy współczynnik korelacji między poszczególnymi kolorami oryginalnego obrazu
# a kolorami zakodowanego obrazu


wspolczynnik_korelacj <- function(C_w, T_C) {
  
  licznik <- sum((C_w - mean(C_w))*(T_C - mean(T_C)))   
  mianownik <- sqrt(sum(((C_w - mean(C_w)))^2))*sqrt(sum(((T_C - mean(T_C)))^2))
  return(licznik/mianownik)
}

corr_R <- wspolczynnik_korelacj(R_w, T_R)
corr_G <- wspolczynnik_korelacj(G_w, T_G)
corr_B <- wspolczynnik_korelacj(B_w, T_B)

# 12. Sprawdzamy wrażliwość kodowania względem doboru klucza.
# przyjmujemy x0 = 0.49, r = 3.79

D_R_zmienione <- dekodowanie(T_R, x0 = 0.49, r = 3.79)
D_G_zmienione <- dekodowanie(T_G, x0 = 0.49, r = 3.79)
D_B_zmienione <- dekodowanie(T_B, x0 = 0.49, r = 3.79)

D_zmienione <- array(c(D_R_zmienione, D_G_zmienione, D_B_zmienione), dim=dim(img))

imageShow(D_zmienione)









