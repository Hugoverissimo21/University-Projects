#2
library(palmerpenguins)
penguins


#21
# selecionar as colunas que queremos
peng_clean <- penguins[, c("species", "sex", "body_mass_g")]

# remover linhas com NA
peng_clean <- na.omit(peng_clean)

# verificar a classe das colunas
cat(class(peng_clean$species), "&", class(peng_clean$sex))

summary(peng_clean)
peng_clean


#22
attach(peng_clean)

# para podermos analisar o fator residuos
npaov <- aov(formula = body_mass_g ~ species * sex, data = peng_clean)

#221
# mudar o layout gráfico para o tipo i,j
par(mfrow = c(1, 2))

# gráfico qq
qqnorm(residuals(npaov))
qqline(residuals(npaov))

# histograma
hist(residuals(npaov), probability = TRUE, ylab = "probability")
xfit <- seq(min(residuals(npaov)), max(residuals(npaov)), length=100)
yfit_residuals <- dnorm(xfit, mean=0, sd = sqrt(var(residuals(npaov))))
lines(xfit, yfit_residuals, col = "black", lwd = 2)


# ver a média pela caixa de bigodes, acho que não vale a pena
#plot(species, residuals(npaov), ylab = "residuals", xlab = "species", cex.axis = 0.7)
#plot(sex, residuals(npaov), ylab = "residuals", xlab = "sex")
# !!!!!!!!!!!!!!!!!!!!!!!!!!!

# voltar ao layout gráfico normal
par(mfrow = c(1, 1))

#norm dos res
shapiro.test(residuals(npaov))

#norm das obs
aggregate(body_mass_g ~ species * sex, data = peng_clean,
          function(x) shapiro.test(x)$p.value)


#222
# independencia dos residuos
#enunciado diz: Assume-se que os erros são independentes.


#FAZER O COISO CHEIO DE PONTOS A TOA NO QUADRADO

# Plot residuals against the predicted values or group levels
plot(residuals(npaov), ylab = "Residuals", xlab = "Observation Index",
     main = "Residuals Plot for ANOVA")
# Add a horizontal line at y=0 for reference
abline(h = 0, col = "red", lty = 1)  

# independencia das observacoes
# H0: não há diferença ou não há associação entre variáveis aka independentes
chisq.test(table(species, sex))

#223
#dispersao entre valores ajustados e residuos
plot(fitted(npaov), residuals(npaov), col = "darkgray", pch = 10)
# linha que melhor se ajusta aos padroes dos residuos
lines(lowess(residuals(npaov) ~ fitted(npaov)), col = "black")
# MUDAR A COR para tons de cinzento preto idk, transparencias das bolas???

# aplicar transformacao logaritmica para tornar padroes mais evidentes
log_resid <- log1p(abs(residuals(npaov)))
plot(fitted(npaov), log_resid, col = "darkgray", pch = 10)
lines(lowess(log_resid ~ fitted(npaov)), col = "black")

detach(peng_clean)
# introduzir coluna nova do tipo "species.sex" (ex.: Gentoo.female)
peng_clean$bart <- interaction(peng_clean$species, peng_clean$sex)
attach(peng_clean)

summary(bart)
# repara que (e.g., nij ≥ 5 em todas as células), e amostra grande, ent Bartlett
# não sei se é assim
# !!!!!!!!!!!!!!!!!!!!!!!!

bartlett.test(body_mass_g ~ bart, data = peng_clean)

detach(peng_clean)


#23
# fazer o grafico de interacao
attach(peng_clean) 
with(peng_clean, interaction.plot(species, sex, body_mass_g,
                                  type = "b", pch = 19, fixed = T,
                                  xlab = "species", ylab = "media body_mass_g"))

with(peng_clean, interaction.plot(sex, species, body_mass_g,
                                  type = "b", pch = 19, fixed = T,
                                  xlab = "sex", ylab = "media body_mass_g"))
# cruzar e n cruzar ns q, interage ou pode interagir bla bla ns
# N É CRUZAR, É SER OU N PARALELO
detach(peng_clean)


#24
summary(npaov)

# comp multiplas species
TukeyHSD(npaov, "species")
plot(TukeyHSD(npaov, "species"), yaxt = "n")
# mudar legendas pq n cabe
custom_labels <- c("Gen-Chi", "Gen-Ade", "Chi-Ade")
axis(2, at = 1:length(custom_labels),
     labels = custom_labels, las = 1, cex.axis = 0.8)

# comp multiplas sex
TukeyHSD(npaov, "sex")
plot(TukeyHSD(npaov, "sex"))



par(mfrow = c(1, 2))

par(cex.lab = 1.3, cex.axis = 1.3, cex=0.6)
#VER MELHOR ISTO DOS TAMANHOS DA LETRA IG

with(peng_clean, interaction.plot(species, sex, body_mass_g,
                                  type = "b", pch = 19, fixed = T,
                                  xlab = "species", ylab = "media body_mass_g", legend = FALSE))
legend("bottomright", legend = c("female", "male"),
       title = "sex", col = 1, lty = c(2,1), pch = 19)

with(peng_clean, interaction.plot(sex, species, body_mass_g,
                                  type = "b", pch = 19, fixed = T,
                                  xlab = "sex", ylab = "media body_mass_g", legend = FALSE))
legend("bottomright", legend = c("Adelie", "Chinstrap", "Gentoo"),
       title = "sex", col = 1, lty = c(3,2,1), pch = 19)


par(mfrow = c(1, 1))


#

#help("interaction.plot")
############################## v6



attach(peng_clean) 


aggregate(body_mass_g ~ bart, data = peng_clean, FUN = mean)

aggregate(body_mass_g ~ species + sex, data = peng_clean, FUN = mean)



par(mfrow = c(1, 2))
# mudar tamanho da fonte (lab, axis, tudo) em %
par(cex.lab = 1.3, cex.axis = 1.3, cex=0.6)

# fazer o grafico de interacao
with(peng_clean, interaction.plot(species, sex, body_mass_g,
                                  type = "b", pch = 19, fixed = T,
                                  xlab = "species", ylab = "media body_mass_g", legend = FALSE))
legend("bottomright", legend = c("female", "male"),
       title = "species", lty = c(2,1), pch = 19)

with(peng_clean, interaction.plot(sex, species, body_mass_g,
                                  type = "b", pch = 19, fixed = T,
                                  xlab = "sex", ylab = "media body_mass_g", legend = FALSE))
legend("bottomright", legend = c("Adelie", "Chinstrap", "Gentoo"),
       title = "sex", lty = c(3,2,1), pch = 19)

# cruzar e n cruzar ns q, interage ou pode interagir bla bla ns
# N É CRUZAR, É SER OU N PARALELO



plot(species, body_mass_g, ylab = "residuals", xlab = "species", cex.axis = 0.7)

plot(sex, body_mass_g, ylab = "residuals", xlab = "species", cex.axis = 0.7)
par(mfrow = c(1, 1))
detach(peng_clean)


################ v7

# comp multiplas species
TukeyHSD(npaov, "species")
plot(TukeyHSD(npaov, "species"), yaxt = "n")
# mudar legendas pq n cabe
custom_labels <- c("Gen-Chi", "Gen-Ade", "Chi-Ade")
axis(2, at = 1:length(custom_labels),
     labels = custom_labels, las = 1, cex.axis = 0.8)

TukeyHSD(npaov, "species:sex")
plot(TukeyHSD(npaov, "species:sex"))#, yaxt = "n")

# organizar grafico melhor com coisas e ns q
# legendas n aparecem de lado

par(cex.lab = 1.3, cex.axis = 1.3, cex=0.6)
TukeyHSD(npaov, "species:sex")
par(mar = c(5, 10, 4, 3))
plot(TukeyHSD(npaov, "species:sex"), cex.axis = 0.8, las=1)





