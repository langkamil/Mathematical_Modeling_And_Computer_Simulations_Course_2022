# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Burczyk Aleksandra
# Gmaj Monika
# Jankowska Kinga 
# Langowski Kamil 
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# PROJEKT 2
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ZADANIE 2.1
# Wyprowadzenie niestandardowej metody Mickensa 
# dla modelu SIR
# dS(t) / dt = m * N - (b / N) * S(t) * I(t) - m * S(t)
# dI(t) / dt = (b / N) * S(t) * I(t) - g * I(t) - m * I(t)
# R(t) = N - S(t) - I(t)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# implementacja metody mickensa 
# wyznaczone 
# S_{n+1}  =  S_n1
# I_{n+1}  =  I_n1 
# R_{n+1}  =  R_n1
# S_{n} = Sn, I_{n} = In, R_{n} = Rn

S_n1 <- function(h, m, N, Sn, b, In) return((h*m*N+Sn)/(1+h*m+(b/N)*h*In))

I_n1 <- function(b, N, h, S_n1, In, g, m) return(((b/N)*h*S_n1*In + In)/(1+g*h+m*h))

R_n1 <- function(g, h, I_n1, Rn, m) return((g*h*I_n1+Rn)/(1+m*h))


# funkcja dokonujaca symulacji metody numerycznej mickensa 
simulation_mickens <- function(S_0, I_0, R_0, h, b,g, m, TT, N) {
  n = TT / h
  
  S={}
  S[1] <- S_0
  
  I={}
  I[1] <- I_0
  
  R={}
  R[1] <- R_0
  
  for(i in 2:n) {
    S_n_1 = S_n1(h, m, N, S[i-1], b, I[i-1])
    S[i] = S_n_1
    
    I_n_1 = I_n1(b, N, h, S[i], I[i-1], g, m)
    I[i] = I_n_1
    
    R_n_1 = R_n1(g, h, I[i], R[i-1], m)
    R[i] = R_n_1
  }
  return(list(S,I,R))
}


# funkcja wyznaczajaca wspolczynnik R0
# zgodnie ze wzorem 
# R0 = b / (m + g)
R_0 <- function(b, g, m) return( b / (g + m))


# funkcja wyznaczajaca punkt rownowagi modelu 
pkt_S <- function(N, R_0) return (N/R_0)

pkt_I <- function(N, m, b, R_0) return ((N * m * (R_0 - 1)) / (b))

pkt_R <- function(N, g, b, R_0) return ((N * g * (R_0 - 1)) / (b))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ZADANIE 2.2 a
# Wyprobowujemy 4 zestawy parametrow, zaczynajac od a
# a) (S0, I0, b, g, m) =  (N - 1, 1, 0.5, 0.8, 0.4),
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# - wyznaczymy punkt rownowagi oraz wartosc R0

#   wyznaczamy wspolczynnik R_0 
#   parametry: b = 0.5, g = 0.8, m = 0.4
R0_a = R_0(b = 0.5, g = 0.8, m= 0.4)
R0_a

#   wyznaczamy punkt rownowagi 
#   parametry b = 0.5, g = 0.8, m = 0.4, N = 10000
rownowagaSIR_a <- c(pkt_S(N = 10000, R0_a), pkt_I(N = 10000, b = 0.5, m= 0.4, R0_a), pkt_R(N = 10000, b = 0.5, g = 0.8, R0_a))
rownowagaSIR_a

# - przedstawimy na jednym wykresie rozwiazanie numeryczne modelu (S, I, R)
#   rozwiazanie numeryczne modelu SIR
x_a = simulation_mickens(9999, 1, 0, h = 0.1, b = 0.5, g = 0.8, m = 0.4, TT = 1000, N = 10000)
S_a = x_a[[1]]
I_a = x_a[[2]]
R_a = x_a[[3]]
df_a = data.frame(S_a, I_a, R_a)

plot(c(1:length(S_a)), S_a, main = "wykres S w czasie", type = "l", col = "skyblue3", xlab = "czas ", ylab = "liczba osob")
legend("right", c("susceptibles"), col = c("skyblue3"), lty = 1, bty = "n")

plot(c(1:length(I_a)), I_a, main = "wykres I w czasie", type = "l", col = "palevioletred3", xlab = "czas ", ylab = "liczba osob")
legend("right", c("infectious"), col = c("palevioletred3"), lty = 1, bty = "n")

plot(c(1:length(R_a)), R_a, main = "wykres R w czasie", type = "l", col = "yellowgreen", xlab = "czas ", ylab = "liczba osob")
legend("right", c("recovered"), col = c("yellowgreen"), lty = 1, bty = "n")


plot(c(1:length(R_a)), R_a, main = "wykres SIR", type = "l", col = "yellowgreen", xlab = "czas ", ylab = "liczba osob")
lines(c(1:length(I_a)), I_a, type = "l", col = "palevioletred3")
lines(c(1:length(S_a)), S_a, type = "l", col = "skyblue3")
legend("right", c("susceptibles","infectious", "recovered"), col = c("skyblue3","palevioletred3", "yellowgreen"), lty = 1, bty = "n")



# - wyznaczymy portrety fazowe dla par I ~ S, S ~ R, I ~ R 
#   wraz z zaznaczonymi punktami poczatkowymi oraz punktami rownowagi 

# I ~ S
plot(I_a, S_a, main = " I ~ S ", type = "l")
points(rownowagaSIR_a[2], rownowagaSIR_a[1], pch = 16, col = "darkorchid1")
points(I_a[1], y = S_a[1], pch = 16, col = "coral2")
legend("topright", c("pkt rownowagi", "pkt poczatkowy"),
       col = c("darkorchid1", "coral2"),cex=0.5, pch=c(16,16))

# S ~ R
plot(S_a, R_a, main = " S ~ R ", type = "l")
points(rownowagaSIR_a[1], rownowagaSIR_a[3], pch = 16, col = "darkorchid1")
points(S_a[1], y = R_a[1], pch = 16, col = "coral2")
legend("topright", c("pkt rownowagi", "pkt poczatkowy"),
       col = c("darkorchid1", "coral2"),cex=0.5, pch=c(16,16))


# I ~ R
plot(I_a, R_a, main = " I ~ R ", type = "l")
points(rownowagaSIR_a[2], rownowagaSIR_a[3], pch = 16, col = "darkorchid1")
points(I_a[1], y = R_a[1], pch = 16, col = "coral2")
legend("topright", c("pkt rownowagi", "pkt poczatkowy"),
       col = c("darkorchid1", "coral2"),cex=0.5, pch=c(16,16))


# - wyznaczymy trojwymiarowy portret fazowy 
#   S ~ I ~ R
# install.packages("plotly")
library("plotly")
plot_ly(df_a, x=df_a$S_a, y=~df_a$I_a, z=~df_a$R_a, main = " S ~ I ~ R ", mode = 'lines', opacity = 1)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ZADANIE 2.2 b
# b) (S0, I0, b, g, m) =  (N - 1, 1, 0.2, 0.04, 0.01)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# - wyznaczymy punkt rownowagi oraz wartosc R0

#   wyznaczamy wspolczynnik R_0 
#   parametry: b = 0.2, g = 0.04, m = 0.01
R0_b = R_0(b = 0.2, g = 0.04, m = 0.01)
R0_b

#   wyznaczamy punkt rownowagi 
#   parametry b = 0.2, g = 0.04, m = 0.01, N = 10000
rownowagaSIR_b <- c(pkt_S(N = 10000, R0_b), pkt_I(N = 10000, b = 0.2, m = 0.01, R0_b), 
                    pkt_R(N = 10000, b = 0.2, g = 0.04, R0_b))
rownowagaSIR_b

# - przedstawimy na jednym wykresie rozwiazanie numeryczne modelu (S, I, R)
#   rozwiazanie numeryczne modelu SIR
x_b = simulation_mickens(9999, 1, 0, h = 0.1, b = 0.2, g = 0.04, m = 0.01, TT = 1000, N = 10000)
S_b = x_b[[1]]
I_b = x_b[[2]]
R_b = x_b[[3]]
df_b = data.frame(S_b, I_b, R_b)

plot(c(1:length(S_b)), S_b, main = "wykres S w czasie", type = "l", col = "skyblue3", xlab = "czas ", ylab = "liczba osob")
legend("right", c("susceptibles"), col = c("skyblue3"), lty = 1, bty = "n")

plot(c(1:length(I_b)), I_b, main = "wykres I w czasie", type = "l", col = "palevioletred3", xlab = "czas ", ylab = "liczba osob")
legend("right", c("infectious"), col = c("palevioletred3"), lty = 1, bty = "n")

plot(c(1:length(R_b)), R_b, main = "wykres R w czasie", type = "l", col = "yellowgreen", xlab = "czas ", ylab = "liczba osob")
legend("right", c("recovered"), col = c("yellowgreen"), lty = 1, bty = "n")

plot(c(1:length(S_b)), S_b, main = "wykres S I R", type = "l", col = "skyblue3", xlab = "czas ", ylab = "liczba osob")
lines(c(1:length(I_b)), I_b, type = "l", col = "palevioletred3")
lines(c(1:length(R_b)), R_b, type = "l", col = "yellowgreen")
legend("topright", c("susceptibles", "infectious", "recovered"), col = c("skyblue3", "palevioletred3", "yellowgreen"), lty = 1, bty = "n",cex=0.5)

# - wyznaczymy portrety fazowe dla par I ~ S, S ~ R, I ~ R 
#   wraz z zaznaczonymi punktami poczatkowymi oraz punktami rownowagi 

# I ~ S
plot(I_b, S_b, main = " I ~ S ", type = "l")
points(rownowagaSIR_b[2], rownowagaSIR_b[1], pch = 16, col = "darkorchid1")
points(I_b[1], y = S_b[1], pch = 16, col = "coral2")
legend("topright", c("pkt rownowagi", "pkt poczatkowy"),
       col = c("darkorchid1", "coral2"),cex=0.5, pch=c(16,16))

# S ~ R
plot(S_b, R_b, main = " S ~ R ", type = "l")
points(rownowagaSIR_b[1], rownowagaSIR_b[3], pch = 16, col = "darkorchid1")
points(S_b[1], y = R_b[1], pch = 16, col = "coral2")
legend("topright", c("pkt rownowagi", "pkt poczatkowy"),
       col = c("darkorchid1", "coral2"),cex=0.5, pch=c(16,16))


# I ~ R
plot(I_b, R_b, main = " I ~ R ", type = "l")
points(rownowagaSIR_b[2], rownowagaSIR_b[3], pch = 16, col = "darkorchid1")
points(I_b[1], y = R_b[1], pch = 16, col = "coral2")
legend("topright", c("pkt rownowagi", "pkt poczatkowy"),
       col = c("darkorchid1", "coral2"),cex=0.5, pch=c(16,16))


# - wyznaczymy trojwymiarowy portret fazowy 
#   S ~ I ~ R
# install.packages("plotly")
library("plotly")
plot_ly(df_b, x=df_b$S_b, y=~df_b$I_b, z=~df_b$R_b, main = " S ~ I ~ R ", mode = 'lines', opacity = 1)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ZADANIE 2.2 c
# c) (S0, I0, b, g, m) =  (N - 1, 1, 0.8, 0.4, 0.01)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# - wyznaczymy punkt rpwnowagi oraz wartosc R0

#   wyznaczamy wspolczynnik R_0 
#   parametry: b = 0.8, g = 0.4, m = 0.01
R0_c = R_0(b = 0.8, g = 0.4, m = 0.01)
R0_c

#   wyznaczamy punkt rownowagi 
#   parametry b = 0.2, g = 0.04, m = 0.01, N = 10000
rownowagaSIR_c <- c(pkt_S(N = 10000, R0_c), pkt_I(N = 10000, b = 0.8, m = 0.01, R0_c), 
                    pkt_R(N = 10000, b = 0.8, g = 0.4, R0_c))
rownowagaSIR_c

# - przedstawimy na jednym wykresie rozwiazanie numeryczne modelu (S, I, R)
#   rozwiazanie numeryczne modelu SIR
x_c = simulation_mickens(9999, 1, 0, h = 0.1, b = 0.8, g = 0.4, m = 0.01, TT = 1000, N = 10000)
S_c = x_c[[1]]
I_c = x_c[[2]]
R_c = x_c[[3]]
df_c = data.frame(S_c, I_c, R_c)

plot(c(1:length(S_c)), S_c, main = "wykres S w czasie", type = "l", col = "skyblue3", xlab = "czas ", ylab = "liczba osob")
legend("right", c("susceptibles"), col = c("skyblue3"), lty = 1, bty = "n")

plot(c(1:length(I_c)), I_c, main = "wykres I w czasie", type = "l", col = "palevioletred3", xlab = "czas ", ylab = "liczba osob")
legend("right", c("infectious"), col = c("palevioletred3"), lty = 1, bty = "n")

plot(c(1:length(R_c)), R_c, main = "wykres R w czasie", type = "l", col = "yellowgreen", xlab = "czas ", ylab = "liczba osob")
legend("right", c("recovered"), col = c("yellowgreen"), lty = 1, bty = "n")

plot(c(1:length(R_c)), R_c, main = "wykres SIR", type = "l", col = "yellowgreen", xlab = "czas ", ylab = "liczba osob")
lines(c(1:length(I_c)), I_c, type = "l", col = "palevioletred3")
lines(c(1:length(S_c)), S_c, type = "l", col = "skyblue3")
legend("right", c("susceptibles","infectious", "recovered"), col = c("skyblue3","palevioletred3", "yellowgreen"), lty = 1, bty = "n", cex=0.5)


# - wyznaczymy portrety fazowe dla par I ~ S, S ~ R, I ~ R 
#   wraz z zaznaczonymi punktami poczatkowymi oraz punktami rownowagi 

# I ~ S
plot(I_c, S_c, main = " I ~ S ", type = "l")
points(rownowagaSIR_c[2], rownowagaSIR_c[1], pch = 16, col = "darkorchid1")
points(I_c[1], y = S_c[1], pch = 16, col = "coral2")
legend("topright", c("pkt rownowagi", "pkt poczatkowy"),
       col = c("darkorchid1", "coral2"),cex=0.5, pch=c(16,16))

# S ~ R
plot(S_c, R_c, main = " S ~ R ", type = "l")
points(rownowagaSIR_c[1], rownowagaSIR_c[3], pch = 16, col = "darkorchid1")
points(S_c[1], y = R_c[1], pch = 16, col = "coral2")
legend("topright", c("pkt rownowagi", "pkt poczatkowy"),
       col = c("darkorchid1", "coral2"),cex=0.5, pch=c(16,16))


# I ~ R
plot(I_c, R_c, main = " I ~ R ", type = "l")
points(rownowagaSIR_c[2], rownowagaSIR_c[3], pch = 16, col = "darkorchid1")
points(I_c[1], y = R_c[1], pch = 16, col = "coral2")
legend("topright", c("pkt rownowagi", "pkt poczatkowy"),
       col = c("darkorchid1", "coral2"),cex=0.5, pch=c(16,16))


# - wyznaczymy trojwymiarowy portret fazowy 
#   S ~ I ~ R
# install.packages("plotly")
library("plotly")
plot_ly(df_c, x=df_c$S_c, y=~df_c$I_c, z=~df_c$R_c, main = " S ~ I ~ R ", mode = 'lines', opacity = 1)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ZADANIE 2.2 d
# d) (S0, I0, b, g, m) =  (N - 1, 1, 0.8, 0.7, 0.01)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# - wyznaczymy punkt rownowagi oraz wartosc R0

#   wyznaczamy wspolczynnik R_0 
#   parametry: b = 0.8, g = 0.7, m = 0.01
R0_d = R_0(b = 0.8, g = 0.7, m = 0.01)
R0_d

#   wyznaczamy punkt rownowagi 
#   parametry b = 0.2, g = 0.04, m = 0.01, N = 10000
rownowagaSIR_d <- c(pkt_S(N = 10000, R0_d), pkt_I(N = 10000, b = 0.8, m = 0.01, R0_d), 
                    pkt_R(N = 10000, b = 0.8, g = 0.7, R0_d))
rownowagaSIR_d

# - przedstawimy na jednym wykresie rozwiazanie numeryczne modelu (S, I, R)
#   rozwiazanie numeryczne modelu SIR
x_d = simulation_mickens(9999, 1, 0, h = 0.1, b = 0.8, g = 0.7, m = 0.01, TT = 1000, N = 10000)
S_d = x_d[[1]]
I_d = x_d[[2]]
R_d = x_d[[3]]
df_d = data.frame(S_d, I_d, R_d)

plot(c(1:length(S_d)), S_d, main = "wykres S w czasie", type = "l", col = "skyblue3", xlab = "czas ", ylab = "liczba osob")
legend("right", c("susceptibles"), col = c("skyblue3"), lty = 1, bty = "n")

plot(c(1:length(I_d)), I_d, main = "wykres I w czasie", type = "l", col = "palevioletred3", xlab = "czas ", ylab = "liczba osob")
legend("right", c("infectious"), col = c("palevioletred3"), lty = 1, bty = "n")

plot(c(1:length(R_d)), R_d, main = "wykres R w czasie", type = "l", col = "yellowgreen", xlab = "czas ", ylab = "liczba osob")
legend("right", c("recovered"), col = c("yellowgreen"), lty = 1, bty = "n")

plot(c(1:length(R_d)), R_d, main = "wykres SIR", type = "l", col = "yellowgreen", xlab = "czas ", ylab = "liczba osob")
lines(c(1:length(I_d)), I_d, type = "l", col = "palevioletred3")
lines(c(1:length(S_d)), S_d, type = "l", col = "skyblue3")
legend("right", c("susceptibles","infectious", "recovered"), col = c("skyblue3","palevioletred3", "yellowgreen"), lty = 1, bty = "n")


# - wyznaczymy portrety fazowe dla par I ~ S, S ~ R, I ~ R 
#   wraz z zaznaczonymi punktami poczatkowymi oraz punktami rownowagi 

# I ~ S
plot(I_d, S_d, main = " I ~ S ", type = "l")
points(rownowagaSIR_d[2], rownowagaSIR_d[1], pch = 16, col = "darkorchid1")
points(I_d[1], y = S_d[1], pch = 16, col = "coral2")
legend("topright", c("pkt rownowagi", "pkt poczatkowy"),
       col = c("darkorchid1", "coral2"),cex=0.5, pch=c(16,16))

# S ~ R
plot(S_d, R_d, main = " S ~ R ", type = "l")
points(rownowagaSIR_d[1], rownowagaSIR_d[3], pch = 16, col = "darkorchid1")
points(S_d[1], y = R_d[1], pch = 16, col = "coral2")
legend("topright", c("pkt rownowagi", "pkt poczatkowy"),
       col = c("darkorchid1", "coral2"),cex=0.5, pch=c(16,16))


# I ~ R
plot(I_d, R_d, main = " I ~ R ", type = "l")
points(rownowagaSIR_d[2], rownowagaSIR_d[3], pch = 16, col = "darkorchid1")
points(I_d[1], y = R_d[1], pch = 16, col = "coral2")
legend("topright", c("pkt rownowagi", "pkt poczatkowy"),
       col = c("darkorchid1", "coral2"),cex=0.5, pch=c(16,16))


# - wyznaczymy trojwymiarowy portret fazowy 
#   S ~ I ~ R
# install.packages("plotly")
library("plotly")
plot_ly(df_d, x=df_d$S_d, y=~df_d$I_d, z=~df_d$R_d, main = " S ~ I ~ R ", mode = 'lines', opacity = 1)




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ZADANIE 3
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# dla podpunktu b) z poprzedniego zadania 
#   parametry: b = 0.2, g = 0.04, m = 0.01

# S_0 = 5000
# I_0 = 5000
x_b3v1 = simulation_mickens(5000, 5000, 0, h = 0.1, b = 0.2, g = 0.04, m = 0.01, TT = 1000, N = 10000)
S_b3v1 = x_b3v1[[1]]
I_b3v1 = x_b3v1[[2]]
R_b3v1 = x_b3v1[[3]]



# S_0 = 7000
# I_0 = 3000
x_b3v2 = simulation_mickens(7000, 3000, 0, h = 0.1, b = 0.2, g = 0.04, m = 0.01, TT = 1000, N = 10000)
S_b3v2 = x_b3v2[[1]]
I_b3v2 = x_b3v2[[2]]
R_b3v1 = x_b3v2[[3]]


# S_0 = 1200
# I_0 = 8800
x_b3v3 = simulation_mickens(1200,8800, 0, h = 0.1, b = 0.2, g = 0.04, m = 0.01, TT = 1000, N = 10000)
S_b3v3 = x_b3v3[[1]]
I_b3v3 = x_b3v3[[2]]
R_b3v3 = x_b3v3[[3]]


# S_0 = 6729
# I_0 = 3271
x_b3v4 = simulation_mickens(6729, 3271, 0, h = 0.1, b = 0.2, g = 0.04, m = 0.01, TT = 1000, N = 10000)
S_b3v4 = x_b3v4[[1]]
I_b3v4 = x_b3v4[[2]]
R_b3v4 = x_b3v4[[3]]


# S_0 = 9900 
# I_0 = 100
x_b3v5 = simulation_mickens(9900, 100, 0, h = 0.1, b = 0.2, g = 0.04, m = 0.01, TT = 1000, N = 10000)
S_b3v5 = x_b3v5[[1]]
I_b3v5 = x_b3v5[[2]]
R_b3v5 = x_b3v5[[3]]


# org, 1, 3, 4, 5
plot(c(1:length(I_b)), I_b, main = "wykres I (zainfekowani) dla roznych warunkow poczatkowych", type = "l", col = "palevioletred3", xlab = "czas ", ylab = "liczba osob")
lines(c(1:length(I_b)),I_b3v1, col="cadetblue")
lines(c(1:length(I_b)),I_b3v3, col="darksalmon")
lines(c(1:length(I_b)),I_b3v4, col="mediumpurple")
lines(c(1:length(I_b)),I_b3v5, col="steelblue")
legend("topright", c("S_0 = 9999, I_0 = 1", 
                     "S_0 = 5000, I_0 = 5000", "S_0 = 1200, I_0 = 8800",  
                     "S_0 = 6729, I_0 = 3271",  "S_0 = 9900, I_0 = 100")
       , col = c("palevioletred3", "cadetblue","darksalmon","mediumpurple",
                 "steelblue" ),
       lty = 1, bty = "n", cex=0.5)


# I S
plot(I_b, S_b, main = " I ~ S ", type = "l", col="palevioletred3", xlim = c(0,10000))
lines(I_b3v1, S_b3v1, col = "cadetblue")
lines(I_b3v3, S_b3v3, col = "darksalmon")
lines(I_b3v4, S_b3v4, col = "mediumpurple")
lines(I_b3v5, S_b3v5, col = "steelblue")
points(rownowagaSIR_b[2], rownowagaSIR_b[1], pch = 16, col = "slategray1")
points(I_b[1], y = S_b[1], pch = 16, col = "palevioletred3")
points(I_b3v1[1], y = S_b3v1[1], pch = 16, col = "cadetblue")
points(I_b3v3[1], y = S_b3v3[1], pch = 16, col = "darksalmon")
points(I_b3v4[1], y = S_b3v4[1], pch = 16, col = "mediumpurple")
points(I_b3v5[1], y = S_b3v5[1], pch = 16, col = "steelblue")
legend("topright", c("pkt rownowagi","S_0 = 9999, I_0 = 1", 
                     "S_0 = 5000, I_0 = 5000", "S_0 = 1200, I_0 = 8800",  
                     "S_0 = 6729, I_0 = 3271",  "S_0 = 9900, I_0 = 100")
       , col = c("slategray1","palevioletred3", "cadetblue","darksalmon","mediumpurple",
                 "steelblue" ),
       lty = 1, bty = "n",cex=0.5, pch=c(16,16), inset = -0.12)


# I R
plot(I_b, R_b, main = " I ~ R ", type = "l", col="palevioletred3", xlim = c(0,10000))
lines(I_b3v1, R_b3v1, col = "cadetblue")
lines(I_b3v3, R_b3v3, col = "darksalmon")
lines(I_b3v4, R_b3v4, col = "mediumpurple")
lines(I_b3v5, R_b3v5, col = "steelblue")
points(rownowagaSIR_b[2], rownowagaSIR_b[3], pch = 16, col = "slategray1")
points(I_b[1], y = R_b[1], pch = 16, col = "palevioletred3")
points(I_b3v1[1], y = R_b3v1[1], pch = 16, col = "cadetblue")
points(I_b3v3[1], y = R_b3v3[1], pch = 16, col = "darksalmon")
points(I_b3v4[1], y = R_b3v4[1], pch = 16, col = "mediumpurple")
points(I_b3v5[1], y = R_b3v5[1], pch = 16, col = "steelblue")
legend("topright", c("pkt rownowagi","S_0 = 9999, I_0 = 1", 
                     "S_0 = 5000, I_0 = 5000", "S_0 = 1200, I_0 = 8800",  
                     "S_0 = 6729, I_0 = 3271",  "S_0 = 9900, I_0 = 100")
       , col = c("slategray1","palevioletred3", "cadetblue","darksalmon","mediumpurple",
                 "steelblue" ),
       lty = 1, bty = "n",cex=0.5, pch=c(16,16), inset = -0.12)







# dla podpunktu c) z poprzedniego zadania 
#   parametry: b = 0.8, g = 0.4, m = 0.01

# S_0 = 5000
# I_0 = 5000
x_c3v1 = simulation_mickens(5000, 5000, 0, h = 0.1, b = 0.8, g = 0.4, m = 0.01, TT = 1000, N = 10000)
S_c3v1 = x_c3v1[[1]]
I_c3v1 = x_c3v1[[2]]
R_c3v1 = x_c3v1[[3]]



# S_0 = 7000
# I_0 = 3000
x_c3v2 = simulation_mickens(7000, 3000, 0, h = 0.1, b = 0.8, g = 0.4, m = 0.01, TT = 1000, N = 10000)
S_c3v2 = x_c3v2[[1]]
I_c3v2 = x_c3v2[[2]]
R_c3v1 = x_c3v2[[3]]


# S_0 = 1200
# I_0 = 8800
x_c3v3 = simulation_mickens(1200,8800, 0, h = 0.1, b = 0.8, g = 0.4, m = 0.01, TT = 1000, N = 10000)
S_c3v3 = x_c3v3[[1]]
I_c3v3 = x_c3v3[[2]]
R_c3v3 = x_c3v3[[3]]


# S_0 = 6729
# I_0 = 3271
x_c3v4 = simulation_mickens(6729, 3271, 0, h = 0.1, b = 0.8, g = 0.4, m = 0.01, TT = 1000, N = 10000)
S_c3v4 = x_c3v4[[1]]
I_c3v4 = x_c3v4[[2]]
R_c3v4 = x_c3v4[[3]]


# S_0 = 9900 
# I_0 = 100
x_c3v5 = simulation_mickens(9900, 100, 0, h = 0.1, b = 0.8, g = 0.4, m = 0.01, TT = 1000, N = 10000)
S_c3v5 = x_c3v5[[1]]
I_c3v5 = x_c3v5[[2]]
R_c3v5 = x_c3v5[[3]]


# org, 1, 3, 4, 5
plot(c(1:length(I_c)), I_c, main = "wykres I (zainfekowani) dla roznych warunkow poczatkowych", type = "l", col = "palevioletred3", xlab = "czas ", ylab = "liczba osob")
lines(c(1:length(I_c)),I_c3v1, col="cadetblue")
lines(c(1:length(I_c)),I_c3v3, col="darksalmon")
lines(c(1:length(I_c)),I_c3v4, col="mediumpurple")
lines(c(1:length(I_c)),I_c3v5, col="steelblue")
legend("topright", c("S_0 = 9999, I_0 = 1", 
                     "S_0 = 5000, I_0 = 5000", "S_0 = 1200, I_0 = 8800",  
                     "S_0 = 6729, I_0 = 3271",  "S_0 = 9900, I_0 = 100")
       , col = c("palevioletred3", "cadetblue","darksalmon","mediumpurple",
                 "steelblue" ),
       lty = 1, bty = "n", cex=0.5)


# I S
plot(I_c, S_c, main = " I ~ S ", type = "l", col="palevioletred3", xlim = c(0,10000), ylim = c(0,10000))
lines(I_c3v1, S_c3v1, col = "cadetblue")
lines(I_c3v3, S_c3v3, col = "darksalmon")
lines(I_c3v4, S_c3v4, col = "mediumpurple")
lines(I_c3v5, S_c3v5, col = "steelblue")
points(rownowagaSIR_c[2], rownowagaSIR_c[1], pch = 16, col = "slategray1")
points(I_c[1], y = S_c[1], pch = 16, col = "palevioletred3")
points(I_c3v1[1], y = S_c3v1[1], pch = 16, col = "cadetblue")
points(I_c3v3[1], y = S_c3v3[1], pch = 16, col = "darksalmon")
points(I_c3v4[1], y = S_c3v4[1], pch = 16, col = "mediumpurple")
points(I_c3v5[1], y = S_c3v5[1], pch = 16, col = "steelblue")
legend("topright", c("pkt rownowagi","S_0 = 9999, I_0 = 1", 
                     "S_0 = 5000, I_0 = 5000", "S_0 = 1200, I_0 = 8800",  
                     "S_0 = 6729, I_0 = 3271",  "S_0 = 9900, I_0 100")
       , col = c("slategray1","palevioletred3", "cadetblue","darksalmon","mediumpurple",
                 "steelblue" ),
       lty = 1, bty = "n",cex=0.5, pch=c(16,16), inset = -0.12)


# I R
plot(I_c, R_c, main = " I ~ R ", type = "l", col="palevioletred3", xlim = c(0,10000), ylim = c(0,10000))
lines(I_c3v1, R_c3v1, col = "cadetblue")
lines(I_c3v3, R_c3v3, col = "darksalmon")
lines(I_c3v4, R_c3v4, col = "mediumpurple")
lines(I_c3v5, R_c3v5, col = "steelblue")
points(rownowagaSIR_c[2], rownowagaSIR_c[3], pch = 16, col = "slategray1")
points(I_c[1], y = R_c[1], pch = 16, col = "palevioletred3")
points(I_c3v1[1], y = R_c3v1[1], pch = 16, col = "cadetblue")
points(I_c3v3[1], y = R_c3v3[1], pch = 16, col = "darksalmon")
points(I_c3v4[1], y = R_c3v4[1], pch = 16, col = "mediumpurple")
points(I_c3v5[1], y = R_c3v5[1], pch = 16, col = "steelblue")
legend("topright", c("pkt rownowagi","S_0 = 9999, I_0 = 1", 
                     "S_0 = 5000, I_0 = 5000", "S_0 = 1200, I_0 = 8800",  
                     "S_0 = 6729, I_0 = 3271",  "S_0 = 9900, I_0 100")
       , col = c("slategray1","palevioletred3", "cadetblue","darksalmon","mediumpurple",
                 "steelblue" ),
       lty = 1, bty = "n",cex=0.5, pch=c(16,16), inset = -0.12)

# dla podpunktu d) z poprzedniego zadania 
#   parametry: b = 0.8, g = 0.7, m = 0.01

# S_0 = 5000
# I_0 = 5000
x_d3v1 = simulation_mickens(5000, 5000, 0, h = 0.1, b = 0.8, g = 0.7, m = 0.01, TT = 1000, N = 10000)
S_d3v1 = x_d3v1[[1]]
I_d3v1 = x_d3v1[[2]]
R_d3v1 = x_d3v1[[3]]



# S_0 = 7000
# I_0 = 3000
x_d3v2 = simulation_mickens(7000, 3000, 0, h = 0.1, b = 0.8, g = 0.7, m = 0.01, TT = 1000, N = 10000)
S_d3v2 = x_d3v2[[1]]
I_d3v2 = x_d3v2[[2]]
R_d3v1 = x_d3v2[[3]]


# S_0 = 1200
# I_0 = 8800
x_d3v3 = simulation_mickens(1200,8800, 0, h = 0.1, b = 0.8, g = 0.7, m = 0.01, TT = 1000, N = 10000)
S_d3v3 = x_d3v3[[1]]
I_d3v3 = x_d3v3[[2]]
R_d3v3 = x_d3v3[[3]]


# S_0 = 6729
# I_0 = 3271
x_d3v4 = simulation_mickens(6729, 3271, 0, h = 0.1, b = 0.8, g = 0.7, m = 0.01, TT = 1000, N = 10000)
S_d3v4 = x_d3v4[[1]]
I_d3v4 = x_d3v4[[2]]
R_d3v4 = x_d3v4[[3]]


# S_0 = 9900 
# I_0 = 100
x_d3v5 = simulation_mickens(9900, 100, 0, h = 0.1, b = 0.8, g = 0.7, m = 0.01, TT = 1000, N = 10000)
S_d3v5 = x_d3v5[[1]]
I_d3v5 = x_d3v5[[2]]
R_d3v5 = x_d3v5[[3]]


# org, 1, 3, 4, 5
plot(c(1:length(I_d)), I_d, main = "wykres I (zainfekowani) dla roznych warunkow poczatkowych", type = "l", col = "palevioletred3", xlab = "czas ", ylab = "liczba osob")
lines(c(1:length(I_d)),I_d3v1, col="cadetblue")
lines(c(1:length(I_d)),I_d3v3, col="darksalmon")
lines(c(1:length(I_d)),I_d3v4, col="mediumpurple")
lines(c(1:length(I_d)),I_d3v5, col="steelblue")
legend("topright", c("S_0 = 9999, I_0 = 1", 
                     "S_0 = 5000, I_0 = 5000", "S_0 = 1200, I_0 = 8800",  
                     "S_0 = 6729, I_0 = 3271",  "S_0 = 9900, I_0 100")
       , col = c("palevioletred3", "cadetblue","darksalmon","mediumpurple",
                 "steelblue" ),
       lty = 1, bty = "n",cex=0.5, pch=c(16,16), inset = -0.12)


# I S
plot(I_d, S_d, main = " I ~ S ", type = "l", col="palevioletred3", xlim = c(0,2000), ylim= c(5000,10000))
lines(I_d3v1, S_d3v1, col = "cadetblue")
lines(I_d3v3, S_d3v3, col = "darksalmon")
lines(I_d3v4, S_d3v4, col = "mediumpurple")
lines(I_d3v5, S_d3v5, col = "steelblue")
points(rownowagaSIR_d[2], rownowagaSIR_d[1], pch = 16, col = "slategray1")
points(I_d[1], y = S_d[1], pch = 16, col = "palevioletred3")
points(I_d3v1[1], y = S_d3v1[1], pch = 16, col = "cadetblue")
points(I_d3v3[1], y = S_d3v3[1], pch = 16, col = "darksalmon")
points(I_d3v4[1], y = S_d3v4[1], pch = 16, col = "mediumpurple")
points(I_d3v5[1], y = S_d3v5[1], pch = 16, col = "steelblue")
legend("topright", c("pkt rownowagi","S_0 = 9999, I_0 = 1", 
                     "S_0 = 5000, I_0 = 5000", "S_0 = 1200, I_0 = 8800",  
                     "S_0 = 6729, I_0 = 3271",  "S_0 = 9900, I_0 100")
       , col = c("slategray1","palevioletred3", "cadetblue","darksalmon","mediumpurple",
                 "steelblue" ),
       lty = 1, bty = "n",cex=0.5, pch=c(16,16), inset = -0.12)


# I R
plot(I_d, R_d, main = " I ~ R ", type = "l", col="palevioletred3", xlim = c(0,10000), ylim=c(0,10000) )
lines(I_d3v1, R_d3v1, col = "cadetblue")
lines(I_d3v3, R_d3v3, col = "darksalmon")
lines(I_d3v4, R_d3v4, col = "mediumpurple")
lines(I_d3v5, R_d3v5, col = "steelblue")
points(rownowagaSIR_d[2], rownowagaSIR_d[3], pch = 16, col = "slategray1")
points(I_d[1], y = R_d[1], pch = 16, col = "palevioletred3")
points(I_d3v1[1], y = R_d3v1[1], pch = 16, col = "cadetblue")
points(I_d3v3[1], y = R_d3v3[1], pch = 16, col = "darksalmon")
points(I_d3v4[1], y = R_d3v4[1], pch = 16, col = "mediumpurple")
points(I_d3v5[1], y = R_d3v5[1], pch = 16, col = "steelblue")
legend("topright", c("pkt rownowagi","S_0 = 9999, I_0 = 1", 
                     "S_0 = 5000, I_0 = 5000", "S_0 = 1200, I_0 = 8800",  
                     "S_0 = 6729, I_0 = 3271",  "S_0 = 9900, I_0 100")
       , col = c("slategray1","palevioletred3", "cadetblue","darksalmon","mediumpurple",
                 "steelblue" ),
       lty = 1, bty = "n",cex=0.5, pch=c(16,16), inset = -0.12)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ZADANIE 4
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# a przeskalowane 
x_a4 = simulation_mickens(9999, 1, 0, h = 0.1, b = 0.5 , g = 0.8, m = 0.4, TT = 1000, N = 10000)
S_a4 = x_a4[[1]]/10000
I_a4 = x_a4[[2]]/10000
R_a4 = x_a4[[3]]/10000

plot(c(1:length(S_a4)), S_a4, main = "wykres S I R (przeskalowane)", type = "l", col = "skyblue3", xlab = "czas ", ylab = "liczba osob", ylim = c(0,1.2))
lines(c(1:length(I_a4)), I_a4, type = "l", col = "palevioletred3")
lines(c(1:length(R_a4)), R_a4, type = "l", col = "yellowgreen")
legend("topright", c("susceptibles", "infectious", "recovered"), col = c("skyblue3", "palevioletred3", "yellowgreen"), lty = 1, bty = "n",cex=0.5)



# b przeskalowane 
x_b4 = simulation_mickens(9999, 1, 0, h = 0.1, b = 0.2, g = 0.04, m = 0.01, TT = 1000, N = 10000)
S_b4 = x_b4[[1]]/10000
I_b4 = x_b4[[2]]/10000
R_b4 = x_b4[[3]]/10000

plot(c(1:length(S_b4)), S_b4, main = "wykres S I R (przeskalowane)", type = "l", col = "skyblue3", xlab = "czas ", ylab = "liczba osob", ylim = c(0,1.2))
lines(c(1:length(I_b4)), I_b4, type = "l", col = "palevioletred3")
lines(c(1:length(R_b4)), R_b4, type = "l", col = "yellowgreen")
legend("topright", c("susceptibles", "infectious", "recovered"), col = c("skyblue3", "palevioletred3", "yellowgreen"), lty = 1, bty = "n",cex=0.5)


# c przeskalowane 
x_c4 = simulation_mickens(9999, 1, 0, h = 0.1, b = 0.8, g = 0.4, m = 0.01, TT = 1000, N = 10000)
S_c4 = x_c4[[1]]/10000
I_c4 = x_c4[[2]]/10000
R_c4 = x_c4[[3]]/10000

plot(c(1:length(S_c4)), S_c4, main = "wykres S I R (przeskalowane)", type = "l", col = "skyblue3", xlab = "czas ", ylab = "liczba osob", ylim = c(0,1.2))
lines(c(1:length(I_c4)), I_c4, type = "l", col = "palevioletred3")
lines(c(1:length(R_c4)), R_c4, type = "l", col = "yellowgreen")
legend("topright", c("susceptibles", "infectious", "recovered"), col = c("skyblue3", "palevioletred3", "yellowgreen"), lty = 1, bty = "n",cex=0.5)


# d przeskalowane 
x_d4 = simulation_mickens(9999, 1, 0, h = 0.1, b = 0.8, g = 0.7, m = 0.01, TT = 1000, N = 10000)
S_d4 = x_d4[[1]]/10000
I_d4 = x_d4[[2]]/10000
R_d4 = x_d4[[3]]/10000

plot(c(1:length(S_d4)), S_d4, main = "wykres S I R (przeskalowane)", type = "l", col = "skyblue3", xlab = "czas ", ylab = "liczba osob", ylim = c(0,1.2))
lines(c(1:length(I_d4)), I_d4, type = "l", col = "palevioletred3")
lines(c(1:length(R_d4)), R_d4, type = "l", col = "yellowgreen")
legend("topright", c("susceptibles", "infectious", "recovered"), col = c("skyblue3", "palevioletred3", "yellowgreen"), lty = 1, bty = "n",cex=0.5)



# b przeskalowane dla roznych liczebnosci populacji
# liczebnosc = 902
x_b4v1 = simulation_mickens(901, 1, 0, h = 0.1, b = 0.2, g = 0.04, m = 0.01, TT = 1000, N = 902)
S_b4v1 = x_b4v1[[1]]/902
I_b4v1 = x_b4v1[[2]]/902
R_b4v1 = x_b4v1[[3]]/902

plot(c(1:length(S_b4v1)), S_b4v1, main = "wykres S I R (przeskalowane), liczebnosc = 902", type = "l", col = "skyblue3", xlab = "czas ", ylab = "liczba osob", ylim = c(0,1.2))
lines(c(1:length(I_b4v1)), I_b4v1, type = "l", col = "palevioletred3")
lines(c(1:length(R_b4v1)), R_b4v1, type = "l", col = "yellowgreen")
legend("topright", c("susceptibles", "infectious", "recovered"), col = c("skyblue3", "palevioletred3", "yellowgreen"), lty = 1, bty = "n",cex=0.5)

# liczebnosc = 5000
x_b4v2 = simulation_mickens(4999, 1, 0, h = 0.1, b = 0.2, g = 0.04, m = 0.01, TT = 1000, N = 5000)
S_b4v2 = x_b4v2[[1]]/5000
I_b4v2 = x_b4v2[[2]]/5000
R_b4v2 = x_b4v2[[3]]/5000

plot(c(1:length(S_b4v2)), S_b4v2, main = "wykres S I R (przeskalowane), liczebnosc = 5 000", type = "l", col = "skyblue3", xlab = "czas ", ylab = "liczba osob", ylim = c(0,1.2))
lines(c(1:length(I_b4v2)), I_b4v2, type = "l", col = "palevioletred3")
lines(c(1:length(R_b4v2)), R_b4v2, type = "l", col = "yellowgreen")
legend("topright", c("susceptibles", "infectious", "recovered"), col = c("skyblue3", "palevioletred3", "yellowgreen"), lty = 1, bty = "n",cex=0.5)

# liczebnosc = 90 000
x_b4v3 = simulation_mickens(89999, 1, 0, h = 0.1, b = 0.2, g = 0.04, m = 0.01, TT = 1000, N = 90000)
S_b4v3 = x_b4v3[[1]]/90000
I_b4v3 = x_b4v3[[2]]/90000
R_b4v3 = x_b4v3[[3]]/90000

plot(c(1:length(S_b4v3)), S_b4v3, main = "wykres S I R (przeskalowane), liczebnosc = 90 000", type = "l", col = "skyblue3", xlab = "czas ", ylab = "liczba osob", ylim = c(0,1.2))
lines(c(1:length(I_b4v3)), I_b4v3, type = "l", col = "palevioletred3")
lines(c(1:length(R_b4v3)), R_b4v3, type = "l", col = "yellowgreen")
legend("topright", c("susceptibles", "infectious", "recovered"), col = c("skyblue3", "palevioletred3", "yellowgreen"), lty = 1, bty = "n",cex=0.5)


# liczebnosc = 1 000 000
x_b4v4 = simulation_mickens(999999, 1, 0, h = 0.1, b = 0.2, g = 0.04, m = 0.01, TT = 1000, N = 1000000)
S_b4v4 = x_b4v4[[1]]/1000000
I_b4v4 = x_b4v4[[2]]/1000000
R_b4v4 = x_b4v4[[3]]/1000000

plot(c(1:length(S_b4v4)), S_b4v4, main = "wykres S I R (przeskalowane), liczebnosc = 1 000 000", type = "l", col = "skyblue3", xlab = "czas ", ylab = "liczba osob", ylim = c(0,1.2))
lines(c(1:length(I_b4v4)), I_b4v4, type = "l", col = "palevioletred3")
lines(c(1:length(R_b4v4)), R_b4v4, type = "l", col = "yellowgreen")
legend("topright", c("susceptibles", "infectious", "recovered"), col = c("skyblue3", "palevioletred3", "yellowgreen"), lty = 1, bty = "n",cex=0.5)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ZADANIE 5
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# instalujemy pakiet
install.packages("epimdr")
library("epimdr")

#wyswietlamy zbior danych niamey
niamey

# sumujemy wartosci cases_1, cases_2, cases_3
suma <- (niamey['cases_1'] + niamey['cases_2'] +niamey['cases_3'])
suma

install.packages("dplyr")
library("dplyr")
# wyswietlamy numery tygodni (1-31)
tyg <- select(niamey, absweek)
tyg

# tworzymy tabelke zawierajaca numer tygodnia oraz liczbe zachorowan dla trzech dystrktow razem
dane <- cbind(tyg, suma)
dane
plot(dane, xlab='tydzien', ylab='liczba zachorowan', main='Liczba zachorowan w poszczegolnych tygodniach dla 3 dystryktow')



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ZADANIE 6
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# zgodnie ze wzorem 
# R0 = b / (m + g)
m = 1/(15*365)
b = 0.79
g = 0.7
R_0_6 <- b / (m+g)
R_0_6

# ponownie przypominamy metode numeryczna mickensa
S_n1 <- function(h, m, N, Sn, b, In) return((h*m*N+Sn)/(1+h*m+(b/N)*h*In))
I_n1 <- function(b, N, h, S_n1, In, g, m) return(((b/N)*h*S_n1*In + In)/(1+g*h+m*h))
R_n1 <- function(g, h, I_n1, Rn, m) return((g*h*I_n1+Rn)/(1+m*h))

# funkcja dokonujaca symulacji metody numerycznej mickensa 
# tutaj zmieniamy R_0, niech bedzie liczone z paramterow b/(m+g), wczesniej bylo argumentem funkcji
simulation_mickens_2 <- function(S_0, I_0, h, b, g, m, TT, N) {
  n = TT / h
  
  S={}
  S[1] <- S_0
  
  I={}
  I[1] <- I_0
  
  R={}
  R[1] <- b/(m+g)
  
  for(i in 2:n) {
    S_n_1 = S_n1(h, m, N, S[i-1], b, I[i-1])
    S[i] = S_n_1
    
    I_n_1 = I_n1(b, N, h, S[i], I[i-1], g, m)
    I[i] = I_n_1
    
    R_n_1 = R_n1(g, h, I[i], R[i-1], m)
    R[i] = R_n_1
  }
  return(list(S,I,R))
}

# podajemy parametry i wartosci wejsciowe
b = 0.79
g = 0.7
m = 1/(15*365)
R_0 = b/(m+g)
N = 769454
I_0 = dane[1,2]
S_0 = N - I_0
h = 5
TT = 1000

symulacja_6 <- simulation_mickens_2(S_0, I_0, h, b, g, m, TT, N) 
symulacja_6

# wybieramy z symulacji druga wartosc, odpowiadaja grupie zarazonych, przeliczamy to na tygodnie
numeryczne <- symulacja_6[[2]][seq(1,length(symulacja_6[[2]]),h)]
numeryczne

plot(numeryczne, xlab='tydzien', ylab='liczba zachorowan', col='red')
points(dane)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ZADANIE 7
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# definiujemy funkcje BNK
par = c(b,g)
BNK <- function(par, h) {
  wyniki <- c()
  I_0 = dane[1,2]
  N = 769454
  S_0 = N - I_0
  m = 1/(15*365)
  TT = 1000
  sym <- simulation_mickens_2(S_0, I_0, h, b, g, m, TT, N)
  num <- sym[[2]][seq(1,length(sym[[2]]),h)]
  for(i in 1:31){
    wyniki[i] <- (dane$cases_1[i] - num[i])^2
  }
  return(sum(wyniki))
}


b = 0.79
g = 0.7
h = 5

BNK(par,h)

