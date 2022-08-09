#MICRO TESTE 3


#Variancia de White . 
#comparar com a matriz de cov ajustada para heterocedasticidade 

##Carregando pacotes necessários
library(car)
library(lmtest)
library(sandwich)

head(Salaries)

##CRIANDO OS OBJETOS: MATRIX X e -----

x0 <- rep(1, 397)
x1 <- matrix(data=Salaries$yrs.since.phd) #anos de estudo
x2 <- matrix(data=Salaries$yrs.service) #idade
x <- matrix(data=cbind(x0, x1, x2), nrow = 397, ncol = 3) #MATRIZ X
y <- matrix(data=Salaries$salary) #salario (Var dependente) MATRIX Y


##cRIANDO A MATRIZ M----

#A matriz M surge das equações normais do MQO após algumas manipulações algébricas.
###Sua fórumula é 

#############m = (I - x(x'x)^-1x').

##Portanto, precisamos criar uma identidade (com traço igual nº de linhas de X).

I <- diag(rep(1, 397))

m = I-(x%*%(solve(t(x)%*%x))%*%t(x))
m       

#RESÍDUOS-----

##Novamente das equações normais, os resíduos são 

############e = my

#Logo, define-se:


e <- m%*%y
head(e)

#TRAÇO DE In e Ik-----

##A variância amostral tem formula

############### s^2 = e'e/(n - k)

###em que n e k são os traços da matriz In e Ik = x(x'x)^-1x.

###Precisamos da variância amostra, pois a variância populacional é desconhecida.

###Portanto, vamos definir esses traços com a função sum do R:

n <- sum(diag(I))
k <- sum(diag(solve(t(x)%*%x)%*%(t(x)%*%x)))


#MATRIZ POSITIVA DEFINIDA -----

#A estratégia aqui foi elevar a matriz-coluna de resíduos ao quadrado e diagona-
#lizá-la para assim obter uma forma funcional da matriz positiva definida da
#variancia de White.

sqr_res <- e^2
mat_diag_res<-diag(as.numeric(sqr_res),397)
View(head(mat_diag_res))

#Variância de White----

#Sua fórmula é:

##########Var_white(B)= (n/(n-k))(X'X)-¹*X'*{sqr_res}*X*(X'X)-¹

Var_white<-((solve(t(x)%*%x))%*%t(x)%*%mat_diag_res%*%x%*%(solve(t(x)%*%x))*(n/(n-k)))
Var_white


#Testar a válidade da matriz encontrada----

##Com o comvando 'vcovHC' do pacotece lmtest verificamos a matriz de covariância
##do estimador de MQO para beta corrigida para heterocedasticidade

vcovHC(lm(y~x1+x2),type="HC1")
