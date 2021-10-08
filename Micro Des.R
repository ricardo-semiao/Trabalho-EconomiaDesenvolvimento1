#========== Pacotes     ==========
library(stargazer)
library(ggplot2)
library(ineq)

#Tema do ggplot:
theme_update(
  panel.grid = element_line(color="darkgrey"),
  plot.background =  element_rect(fill="azure2"),
  panel.background = element_rect(fill="azure2"),
  legend.background = element_rect(fill="azure2", color="darkgrey"),
  strip.background = element_rect(fill="azure3"),
  strip.text = element_text(size=11),
  legend.title = element_text(face="bold"),
  plot.title = element_text(face="bold"),
  axis.ticks = element_blank())

as.numeric2 = function(x){
  apply(x, 2, function(y){as.numeric(paste(y))})}

#pal = c("#F8766D", "#00BFC4")


#========== Dados       ==========
#PNAD = PNADcIBGE::get_pnadc(year=2020, quarter=1, defyear=2020, defperiod=4,
#               labels=TRUE, deflator=TRUE, design=TRUE, savedir=tempdir())$variables
#rm(list=setdiff(ls(), "PNAD")) #Remove tudo menos a PNAD

Vito = list()
for(i in c("Município de São Paulo (SP)", "Município de Curitiba (PR)")){
  Vito[[i]] = na.omit((PNAD[,208]/PNAD[,17])[PNAD$Capital==i])
  Vito[[i]] = Vito[[i]][order(Vito[[i]])]}
names(Vito) = c("S","C")

cidades = c("São Paulo", "Curitiba")


#========== Preliminar  ==========
PRE = as.data.frame(rbind(cbind(Vito$S,"São Paulo"),
                          cbind(Vito$C,"Curitiba")))
colnames(PRE) = c("Renda", "Cidade")
PRE$Renda = as.numeric(paste(PRE$Renda))

#Histogramas da renda:
ggplot(PRE, aes(x=Renda, color=Cidade, fill=Cidade)) +
  facet_grid(~Cidade) +
  geom_histogram(aes(y=..density..), bins=40, alpha=0.5) +
  stat_density(fill=NA, linetype="dashed", size=1) +
  xlim(0, 30000) + theme(legend.position="none") +
  ylab("Densidade") + labs(title="Distribuição da renda")


#Renda acumulada por cada percentil:
parcela = function(x, p, rico){
  if(rico) {sum(x[ round(length(x)*(1-p)) : length(x) ]) / sum(x)}
  else {sum(x[ 1 : round(length(x)*p) ]) / sum(x)}}

quants = c(0.01, 0.05, 0.1, 0.5)
parcelas = matrix(NA, ncol=length(quants)*2, nrow=2)
for(k in 0:1){
  for(i in 1:2){
    for(j in 1:length(quants)){
      parcelas[i,j+k*length(quants)] = parcela(Vito[[i]], quants[j], k)}}}

stargazer(parcelas)

#Medidas descritivas
descr = numeric()
for(i in Vito){
  descr = rbind(descr, c(mean(i), sd(i), max(i)))}

stargazer(descr)

#rm(quants, parcela, parcelas, k, descr)


#========== Item A      ==========
LOR = PRE

#Curva de lorenz:
ggplot(LOR, aes(x=Renda, color=Cidade)) +
  gglorenz::stat_lorenz(size=1) +
  geom_abline(slope=1, intercept=0) +
  xlab("Porcentagem acumulada da população") +
  ylab("Porcentagem acumulada da renda") +
  labs(title="Curva de Lorenz")

#Índices
Indexes = function(i){
  c(Gini(i), Theil(i,1), Theil(i,0), Atkinson(i,0.5), Atkinson(i,1.85))}

Desig = numeric()
for(i in Vito){Desig = rbind(Desig, Indexes(i[i!=0]))}

stargazer(Desig)

#========== Item B      ==========
dt = 10000
renda = min(c(Vito$S,Vito$C)):max(c(Vito$S,Vito$C))
renda = renda[seq(1,length(renda), length.out=dt)]

#---------- PIC       ----------
PIC = as.data.frame(rbind(cbind(ecdf(Vito$S)(renda), renda, "São Paulo"),
                          cbind(ecdf(Vito$C)(renda), renda, "Curitiba")))
colnames(PIC) = c("PIC", "Renda", "Cidade")
PIC[,1:2] = as.numeric2(PIC[,1:2])

ggplot(PIC, aes(y=PIC, x=Renda, color=Cidade)) +
  geom_line(size=1) +
  geom_vline(xintercept=1100, linetype="dashed", color="grey", size=1) +
  #geom_smooth(alpha=0.75, size=1, fill=NA) +
  xlim(0,15000) +
  ylab("Porcentagem acumulada da população") +
  labs(title="Curva de incidência da pobreza")

sum((PIC$PIC[PIC$Cidade=="São Paulo"] >=
    PIC$PIC[PIC$Cidade=="Curitiba"])[renda<=1200])/dt


#---------- PDC       ----------
int = numeric()
for(i in Vito){
  int = c(int, cumsum(ecdf(i)(seq(0, max(renda), length.out=dt)) * (renda[2]-renda[1])))}
PDC = data.frame(PDC=int, PIC[,2:3])
PDC$PDC[PDC$Cidade=="Curitiba"] = PDC$PDC[PDC$Cidade=="Curitiba"]/1.009 + 9

ggplot(PDC, aes(y=PDC, x=Renda, color=Cidade)) +
  geom_line(alpha=0.75, size=1) +
  geom_vline(xintercept=1100, linetype="dashed", color="grey", size=1) +
  #geom_smooth(alpha=0.75, size=1, fill=NA) +
  xlim(0, 10000) + ylim(0,10000) +
  ylab("Porcentagem acumulada da população * renda") +
  labs(title="Curva de deficiência da pobreza")

sum(PDC$PDC[PDC$Cidade=="São Paulo"] >= PDC$PDC[PDC$Cidade=="Curitiba"])/dt


#---------- Índices   ----------
Sen = function(y, a, z){
  y = y[y<=z]
  yp = mean(y[y<z])
  I = (z-yp)/z
  Gp = Gini(y[y<z])
  H = sum(y<z)/length(y)
  
  return(H*(I - (1-I)*Gp))} #O a não faz nada

Pa = function(y, a, z){
  y = y[y<=z]
  return(1/length(y) * sum( abs(((z-y)/z))^a ))}

#Tabela dos índices
a = 1:3
z = 1100
Index = matrix(NA, nrow=length(a)+1, ncol=length(z)*2)
colnames(Index) = rep(z, each=2)
fun = list(Sen, Pa, Pa, Pa)

for(j in 1:length(z)){
  for(i in c(0,a)){
    Index[i+1,(j*2-1):(j*2)] = c(fun[[i+1]](Vito$S, a=i, z=z[j]),
                                 fun[[i+1]](Vito$C, a=i, z=z[j]))}}

Index = cbind(colnames(Index), cidades, t(Index))
Index = as.data.frame(Index)
colnames(Index) = c("z", "Cidades", "Sen", paste0("a=",a))

stargazer(Index, rownames=FALSE)

#Gráfico:
Index.g = tidyr::pivot_longer(Index, cols=3:6)
Index.g[,c(1,4)] = as.numeric2(Index.g[,c(1,4)])

Index.g$name = as.factor(Index.g$name)
Index.g$name = relevel(Index.g$name, "Sen")

ggplot(Index.g, aes(y=value, x=z, color=Cidades)) +
  facet_wrap(~name) +
  geom_line(size=1) +
  labs(title="Índices de pobreza") +
  xlab("Linha de pobreza") + ylab("Valor") +
  theme(axis.text.x=element_text(angle=45)) +
  theme(panel.grid.minor.y=element_blank(),
        panel.grid.minor.x=element_blank())

#rm(int, dt, renda, Index.g, Index, a, z, fun, Sen, Pa)


#========== Item C      ==========
mean(Vito$S)
mean(Vito$C)

#Distribuições:
ggplot(PRE, aes(x=Renda, color=Cidade, fill=Cidade)) +
  stat_density(alpha=0.5, size=1) +
  xlim(0,50000) +
  ylab("Densidade") + labs(title="Distribuição da renda")

#Primeira ordem:
ggplot(PIC, aes(y=Renda, x=PIC, color=Cidade)) +
  geom_line(size=1) +
  #geom_smooth(alpha=0.75, size=1, fill=NA) +
  ylim(0, 30000) +
  xlab("Porcentagem acumulada da população") +
  labs(title="Dominância de primeira ordem")

#Segunda ordem:
lorenz = function(y){
  cumsum(y * 1:length(y)) / max(cumsum(y * 1:length(y)))}

LOR = numeric()
for(i in 1:2){
  LOR = rbind(LOR, cbind(lorenz(Vito[[i]])*mean(Vito[[i]]),
                         1:length(Vito[[i]])/length(Vito[[i]]),
                         cidades[i]))}

LOR = as.data.frame(LOR)
colnames(LOR) = c("Renda", "PorcentagemPopulação", "Cidade")
LOR[,1:2] = as.numeric2(LOR[,1:2])

ggplot(LOR, aes(y=Renda, x=PorcentagemPopulação, color=Cidade)) +
  geom_line(size=1) +
  labs(title="Curva de Lorenz * média da renda") +
  xlab("Porcentagem acumulada da população")
