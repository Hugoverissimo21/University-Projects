############################## FEITO
#2. O conjunto de dados “sat” do package “faraway” foi obtido com o objetivo de estudar
#a relação entre as despesas dos alunos com a educação no ensino público e os resultados
#obtidos no exame SAT. 
#
#install.packages("faraway")
library(faraway)
sat
#
#O conjunto de dados contém 7 variáveis relativas aos resultados de 50 alunos. 

############################## FEITO

# y - SAT (total)
# x - despesas (expend), ordenado (salary), a razão média de alunos por professor (ratio)
#e a percentagem de alunos elegíveis para fazerem o exame (takers).
#
sat_clean <- sat[, c("total", "expend", "salary", "ratio", "takers")]
sat_clean

############################## FEITO

#Ajuste um modelo de regressão linear múltipla aos dados
sat_lm <- lm(total ~ expend + salary + ratio + takers, data = sat_clean) #definir regressaos
sat_lm
summary(sat_lm)
#
#teste a sua significância, alpha = 0,05.
# h0: b1 = b2 = ... = 0 vs h1: existe bi != 0
# f0 = 52.88; p-value < 2.2e-16
# rej h0 para alpha = 0.05,
# i.e., existe evid esta de que o modelo ajustaado é significativo, pelo menos uma
# das variaveis indepen. é signif (bi != 0)
#VI DA AULA, CONFIRMAR
#
# ????
#podemos ver de vez os p-values individuais do tipo da AULA
# h0: b1 = 0 vs h1: b1 != 0
# t0 = 7.532 ; p-value = 0.000653
# rej h0 para alpha = 0.05 , a var. independente x1 é significativa para um modelo
# que inclui a var x2 e a x3 e ...

##### REDISUDOS POR FAZER #########################

# Análise de resíduos para o modelo completo
par(mfrow = c(2, 2))
plot(sat_lm)
#plot(fitted(sat_lm),residuals(sat_lm),ylab="Residuals", xlab="fitted values") # é o primeiro
#abline(0, 0)
#qqnorm(residuals(sat_lm)) # é o segundo
#qqline(residuals(sat_lm))

# e pressupostos do modelo inicial ??

#normalidade
shapiro.test(residuals(sat_lm))

###########??? POR FAZER###################

#Realize todos os testes que achar relevantes para o conjunto de dados em análise.

# ????? ns se estao certos, foi o bot
# Teste de significância global do modelo
anova(sat_lm)

# Testes de significância individual dos coeficientes
summary(sat_lm)$coefficients

############################## FEITO

#Utilize a função stepwise para obter o modelo q ue “melhor” se ajusta ao conjunto de dados “sat”.
sat_lm_step <- step(sat_lm,direction="both")
# menor aic melhor acho eu
# na primeira tabela vemos q fica melhor sem o x4
# na segunda, apos tirar o x4, vemos que estamos no melhor modelo pq todos os aic aumentam
sat_lm_step
summary(sat_lm_step)
# e olha q o Adjusted R-squared:  0.7747  é 77%, aumentou do q x4
# o x4 é mesmo lixo tipo
#
# o y esta relacionado com x1,x2,x3

##############################

#Utilize ainda testes F-parciais para confirmar o resultado obtido pela função stepwise.

# Teste F-parcial
teste_f_parcial <- anova(sat_lm, sat_lm_step)
teste_f_parcial

# Não há evidências suficientes para concluir que as variáveis excluídas são estatisticamente
# importantes para explicar a variação na variável dependente.
# O modelo mais simples é considerado suficiente para explicar os dados,
# pelo menos com base nos critérios utilizados no teste F-parcial.

#POR CONTINUAR ou sao as coisas logo da tabela do T?????

##############################

#Explore da maneira que achar conveniente as funções do R para a regressão linear múltipla.

# PRESSUPOSTOS do modelo final ????

# Identificação de observações influentes
influence.measures(sat_lm_step)

# Matriz de correlação entre as variáveis independentes
cor(sat_clean[, c("salary", "ratio", "takers","expend")])



############## v3

sat_clean2 <- lm(total ~ expend + salary + ratio, data = sat_clean)
#teste_f_parcial <- anova(sat_lm, sat_clean2)
#2.607e-16

sat_clean3 <- lm(total ~ expend + salary  + takers, data = sat_clean)
#teste_f_parcial <- anova(sat_lm, sat_clean3)
#0.2657

sat_clean4 <- lm(total ~ expend + ratio + takers, data = sat_clean)
#teste_f_parcial <- anova(sat_lm, sat_clean4)
#0.4962

sat_clean5 <- lm(total ~ salary + ratio + takers, data = sat_clean)
#teste_f_parcial <- anova(sat_lm, sat_clean5)
#0.6742

sat_clean6 <- lm(total ~ expend + salary, data = sat_clean)
#teste_f_parcial <- anova(sat_lm, sat_clean6)
#1.243e-15

sat_clean7 <- lm(total ~ ratio + takers, data = sat_clean)
#teste_f_parcial <- anova(sat_lm, sat_clean7)
#0.04761

sat_clean8 <- lm(total ~ expend + ratio , data = sat_clean)
#teste_f_parcial <- anova(sat_lm, sat_clean8)
#3.707e-16

sat_clean9 <- lm(total ~  salary + ratio , data= sat_clean)
#teste_f_parcial <- anova(sat_lm, sat_clean9)
#1.491e-15
#teste_f_parcial


min.model <- lm(total ~ 1, data = sat_clean)

step(sat_lm, direction = "both")
step(min.model, direction = "both", scope=formula(sat_lm))
#
step(sat_lm, direction = "backward")
#
step(min.model, direction='forward', scope=formula(sat_lm))

summary(step(min.model, direction='forward', scope=formula(sat_lm)))
summary(sat_lm_step)


max.model <- sat_lm

sat_lm_step_max <- step(max.model, direction = "both")

sat_lm_step_min <- step(min.model, direction = "both", scope = formula(max.model))


#

sat_lm_step_max
sat_lm_step_min

#

summary(sat_lm_step_max)
summary(sat_lm_step_min)



##### v4

#testar h0:b1=0
#> lm2<-lm(formula = y ~ x2, data = lmdata)
#> lm12<-lm(formula = y ~ x1+ x2, data = lmdata)
#> anova(lm2,lm12)


sat_lm

sat0 <- lm(total ~ 1, data = sat_clean)

sat01 <- lm(total ~ expend, data = sat_clean)

summary(sat01)


sat02 <- lm(total ~ salary, data = sat_clean)
summary(sat02)
sat03 <- lm(total ~ ratio, data = sat_clean)
summary(sat03)
sat04 <- lm(total ~ takers, data = sat_clean)
summary(sat04)

anova1 <- anova(sat0,sat01)
anova2 <- anova(sat0,sat02)
anova3 <- anova(sat0,sat03)
anova4 <- anova(sat0,sat04)
anova1


anova1$"Pr(>F)"[2] #p-value
anova2$"Pr(>F)"[2] #p-value
anova3$"Pr(>F)"[2] #p-value
anova4$"Pr(>F)"[2] #p-value


testes_f_parcial <- function(formula_0, variaveis) {
  # modelo inicial
  modelo_0 <- lm(as.formula(formula_0), data = sat_clean)
  
  # fazer anova com a adicao de cada variável
  resultados_anova <- lapply(variaveis, function(x) {
    formula_i <- paste(formula_0, "+", x)
    modelo_i <- lm(as.formula(formula_i), data = sat_clean)
    anova(modelo_0, modelo_i)})
  
  # selecionar os p-values
  p_values <- sapply(resultados_anova, function(resultado) resultado$"Pr(>F)"[2])
  
  # criar dataframe para melhorar a visualizacao
  df <- t(p_values)
  colnames(df) <- variaveis
  rownames(df) <- c("p-values")
  
  return(df)
}

testes_f_parcial("total ~ 1", c("expend", "salary", "ratio", "takers"))

testes_f_parcial("total ~ takers", c("expend", "salary", "ratio"))

testes_f_parcial("total ~ takers + expend", c("salary", "ratio"))



#### v6

summary(sat_lm_step_min)
summary(sat_lm_step_max)

summary(sat_lm_step_min)$adj.r.squared
summary(sat_lm_step_max)$adj.r.squared

sat_lm.step <- sat_lm_step_min



############## v7.0

### v7

#summary(sat_lm.step)
#
summary(sat_lm.step)$call
#
summary(sat_lm.step)$coefficients
#
summary(sat_lm.step)$adj.r.squared
#
summary(sat_lm.step)$fstatistic

# Matriz de correlação entre as variáveis independentes
cor(sat_clean[, c("takers", "expend")])

# Gráfico de resíduos versus fatores
par(mfrow = c(1, 2))
#
# Gráfico para 'expend'
plot(sat_clean$expend, residuals(sat_lm.step), 
     xlab = "expend", ylab = "Residuals",
     main = "Residuals vs. expend")
#
# Gráfico para 'takers'
plot(sat_clean$takers, residuals(sat_lm.step), 
     xlab = "takers", ylab = "Residuals",
     main = "Residuals vs. takers")
#
# Restaurar a configuração padrão
par(mfrow = c(1, 1))

#Um VIF acima de 5 é geralmente considerado um sinal de multicolinearidade significativa,
#e valores acima de 10 indicam multicolinearidade séria.
vif(sat_lm.step)

library(performance, warn.conflicts = FALSE)
check_model(sat_lm.step)
#
#Linearity of the relationships between the dependent and independent variables11
#Independence of the observations
#Normality of the residuals
#Homoscedasticity of the residuals
#No influential points (outliers)
#No multicollinearity
#
#https://statsandr.com/blog/multiple-linear-regression-made-simple/


# se so sao duas variaveis vamos desenhar o grafico !!!!


#library(scatterplot3d)
#library(rgl)

# Assuming 'sat_lm.step' is your linear regression model
# Assuming 'sat_clean' is your data
# Assuming 'total' is your dependent variable and 'takers' and 'expend' are your independent variables

# Make predictions using the model
#predictions <- predict(sat_lm.step, newdata = sat_clean)

# Compute residuals
#residuals <- sat_clean$total - predictions

# Create a 3D scatterplot
#scatterplot3d(
#  x = sat_clean$takers, 
#  y = sat_clean$expend, 
#  z = residuals, 
#  main = "Residuals and Fitted Plane",
#  xlab = "Takers", 
#  ylab = "Expenditure", 
#  zlab = "Residuals"
#)

# Add the fitted plane to the plot
#fit_plane <- sat_lm.step$coefficients[1] + 
#  sat_lm.step$coefficients[2] * sat_clean$takers + 
#  sat_lm.step$coefficients[3] * sat_clean$expend

#plane3d(sat_clean$takers, sat_clean$expend, fit_plane, col = "red", alpha = 0.5)

####
#library(lattice)

# Assuming 'sat_lm.step' is your linear regression model
# Assuming 'sat_clean' is your data
# Assuming 'total' is your dependent variable and 'takers' and 'expend' are your independent variables

# Make predictions using the model
predictions <- predict(sat_lm.step, newdata = sat_clean)

# Compute residuals
residuals <- sat_clean$total - predictions

# Create a 3D scatterplot with regression plane
#cloud(sat_clean$total ~ sat_clean$takers * sat_clean$expend, 
#      groups = NULL, 
#      pch = 16,
#      main = "Residuals and Fitted Plane",
#      xlab = "Takers",
#      ylab = "Expenditure",
#      zlab = "Residuals",
#      panel.3d.cloud = panel.3dscatter,
#      panel.3d.wireframe = function(x, y, z, subscripts, ...) {
#        panel.3dwire(x = x, y = y, z = z, col = "red", ...)
#      },
#      screen = list(z = 35, x = -65),
#      scales = list(arrows = FALSE)
#)




###### v7.2
#summary(sat_lm.step)
#
summary(sat_lm.step)$call
#
summary(sat_lm.step)$coefficients
#
summary(sat_lm.step)$adj.r.squared
#
summary(sat_lm.step)$fstatistic

# Matriz de correlação entre as variáveis independentes
cor(sat_clean[, c("takers", "expend")])

# Gráfico de resíduos versus fatores
par(mfrow = c(1, 2))
#
# Gráfico para 'expend'
plot(sat_clean$expend, residuals(sat_lm.step), 
     xlab = "expend", ylab = "Residuals",
     main = "Residuals vs. expend")
#
# Gráfico para 'takers'
plot(sat_clean$takers, residuals(sat_lm.step), 
     xlab = "takers", ylab = "Residuals",
     main = "Residuals vs. takers")
#
# Restaurar a configuração padrão
par(mfrow = c(1, 1))

#Um VIF acima de 5 é geralmente considerado um sinal de multicolinearidade significativa,
#e valores acima de 10 indicam multicolinearidade séria.
vif(sat_lm.step)

library(performance, warn.conflicts = FALSE)
check_model(sat_lm.step)
#
#Linearity of the relationships between the dependent and independent variables11
#Independence of the observations
#Normality of the residuals
#Homoscedasticity of the residuals
#No influential points (outliers)
#No multicollinearity
#
#https://statsandr.com/blog/multiple-linear-regression-made-simple/

#library(rgl) #nao funciona no Rmd ?








#grafico

# Install and load the scatterplot3d package
library(scatterplot3d)

# Create sample data
x <- c(1, 2, 10,1)
y <- c(1, 2, 10,1)
z <- c(7, 8, 10,1)

# Create a 3D scatter plot
s3d <- scatterplot3d(x, y, z, color="red")#, color = "red", pch = 16, main = "3D Scatter Plot")


# Optionally, you can add more points and lines

# Add more points
x2 <- c(2, 3, 4,0.1)
y2 <- c(5, 6, 7,0.1)
z2 <- c(8, 9, 10,0.1)
s3d$points3d(x2, y2, z2, col="green")#, color = "green", pch = 16, main = "3D Scatter Plot")


# Connect the points with lines
#segments3d(x[3], y[3], z[3], x2[1], y2[1], z2[1], col = "blue", lwd = 2)




#s3d <- scatterplot3d(drill1[1:400], drill1[7:406], drill1[32:431],
#                     color = "red", type = "l", angle = 120, xlab = "drilling torque",
#                     ylab = "drilling torque, lag 6", zlab = "drilling torque, lag 31",
#                     main = "Two deep hole drilling processes")
#s3d$points3d(drill2[1:400], drill2[7:406], drill2[32:431],
#             col = "blue", type = "l")
#legend(s3d$xyz.convert(-400, 1000, 950), col= c("blue", "red"),
#       legend = c("regular process", "chattering process"), lwd = 2,
#       bg = "white")





library(scatterplot3d)


x_min <- min(sat_clean$takers) * 0.9
x_max <- max(sat_clean$takers) * 1.1
y_min <- min(sat_clean$expend) * 0.9
y_max <- max(sat_clean$expend) * 1.1



predictions <- predict(sat_lm.step,
                       newdata = data.frame(takers = c(x_min, x_max), expend = c(y_min, y_max)))


sat_lm.graph <- scatterplot3d(c(x_min,x_max), c(y_min,y_max), predictions,
                     color = "black", type = "l", angle = 30,
                     xlab = "takers", ylab = "expend", zlab = "total",
                     zlim = c(min(sat_clean$total)*0.9, max(sat_clean$total)*1.1),
                     main = "Regressão sat_lm.step")

sat_lm.graph$points3d(sat_clean$takers, sat_clean$expend, sat_clean$total,
             col = "darkgray", type = "p")

summary(sat_lm.step)$fstatistic
summary(sat_lm.step)

summary(sat_lm.step)$fstatistic["Pr(>F)"]



f_statistic <- 106.7
df_between <- 2
df_within <- 47

# Calculate the p-value
pf(f_statistic, df_between, df_within, lower.tail = FALSE)







## v9


# Load required libraries
library(performance, warn.conflicts = FALSE)
library(ggplot2)

# Assuming 'sat_lm.step' is your linear regression model
# Generate the entire set of diagnostic plots
check_model(sat_lm.step)[4]

qqnorm(residuals(sat_lm.step))
qqline(residuals(sat_lm.step))

