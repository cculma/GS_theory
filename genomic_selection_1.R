# load libraries 
library(rrBLUP)
library(e1071)
library(caret)
library(dplyr)
library(ranger)
library(nnet)
library(brnn)
library(arm)
library(monomvn)
library(parallel)
library(doParallel)
library(iterators)
library(foreach)

#to use parallel
cl <- makePSOCKcluster(5)
registerDoParallel(cl)
stopCluster(cl)


setwd("~/OneDrive - Washington State University (email.wsu.edu)/roza_report_2020")
# setwd("/scratch/user/cesar.medinaculma/20201216_115137/GS")

logan <- as.matrix(read.table("NMGS_NGSEP1.2_poly.txt", header = TRUE, row.names = 1))
pheno <- as.matrix(read.csv("BLUPS_YEARS2.csv", header = TRUE, row.names = 1))
dim(logan)
# logan[1:4,1:4]
logan <- logan[,c(1:200)]
pheno <- pheno[,c(1:6)]
head(pheno)
dim(pheno)
col1 <- colnames(pheno)

logan_list <- list()
for (i in 1:ncol(pheno)) {
  j <- as.matrix(pheno[,i][!is.na(pheno[,i])])
  colnames(j) <- "BLUPs"
  logan_list[[length(logan_list)+1]] = j
}

names(logan_list) <- col1

#define training control
train.control <- trainControl(method = 'cv', number = 10, savePredictions = 'final')
set.seed(123)

ranger_list <- list()
svm_list <- list()
bayesglm_list <- list()

for (i in 1:length(logan_list)){
  have.both = intersect(rownames(logan), rownames(logan_list[[i]]))
  logan2 = logan[have.both,]
  pheno1 = logan_list[[i]][have.both,]
  logan3 <- cbind(pheno1, logan2)
  model.ranger <- train(pheno1 ~., data = logan3, method = 'ranger', trControl = train.control)
  model.svm <- train(pheno1 ~., data = logan3, method = 'svmRadialSigma', trControl = train.control)
  model.bayesglm <- train(pheno1 ~., data = logan3, method = 'bayesglm', trControl = train.control)
  
  df.model.ranger <- as.data.frame(print(model.ranger))
  df.model.svm <- as.data.frame(print(model.svm))
  df.model.bayesglm <- as.data.frame(print(model.bayesglm))
  
  model.ranger.all <- df.model.ranger[which.min(df.model.ranger$RMSE),]
  model.svm.all <- df.model.svm[which.min(df.model.svm$RMSE),]
  model.bayesglm.all <- df.model.bayesglm[which.min(df.model.bayesglm$RMSE),]
  
  ranger_list[[length(ranger_list)+1]] = model.ranger.all
  svm_list[[length(svm_list)+1]] = model.svm.all
  bayesglm_list[[length(bayesglm_list)+1]] = model.bayesglm.all
}

names(ranger_list) <- col1
names(svm_list) <- col1
names(bayesglm_list) <- col1

have.both = intersect(rownames(logan), rownames(logan_list[[1]]))
logan2 = logan[have.both,]
pheno1 = logan_list[[5]][have.both,]
dim(logan2)
dim(pheno1)
logan3 <- cbind(pheno1, logan2)
model.ranger <- train(pheno1 ~., data = logan3, method = 'ranger', trControl = train.control)
print(model.ranger)
model.svm <- train(pheno1 ~., data = logan3, method = 'svmRadialSigma', trControl = train.control)
print(model.svm)


df.model.ranger <- as.data.frame(print(model.ranger))
model.ranger.all <- df.model.ranger[which.min(df.model.ranger$RMSE),]
model.ranger.all <- model.ranger.all %>% mutate(pearson = cor(model.ranger[['pred']]$pred, model.ranger[['pred']]$obs),
                                                model = 'ranger',
                                                trait = "September_2019")

(pearson = cor(model.ranger[['pred']]$pred, model.ranger[['pred']]$obs))
length(model.ranger[['pred']]$pred)
length(model.ranger[['pred']]$obs)

obs <- model.ranger[['pred']]$obs
pred <- model.ranger[['pred']]$pred
df <- data.frame(obs, pred, line)
line <- have.both
rownames(df) <- have.both
library(plotly)

fig <- plot_ly(data = df, x = ~obs, y = ~pred) %>% add_markers(alpha = 0.1, ids = ~line)
fig <- plot_ly(data = df, x = ~obs, y = ~pred, text=~line) %>% add_markers(alpha = 0.5)

fig


# model.mlp <- train(pheno1 ~., data = logan3, method = 'mlpWeightDecay', trControl = train.control)
# model.bayesglm <- train(pheno1 ~., data = logan3, method = 'bayesglm', trControl = train.control)
# model.bridge <- train(pheno1 ~., data = logan3, method = 'bridge', trControl = train.control)

save.image(file='ML_GS.RData')


library(plotly)
mtcars %>%
  highlight_key(~cyl) %>%
  plot_ly(
    x = ~wt, y = ~mpg, text = ~cyl, mode = "markers+text", 
    textposition = "top", hoverinfo = "x+y"
  ) %>%
  highlight(on = "plotly_hover", off = "plotly_doubleclick")

df <- read.table("pred_mlp.txt", header = T)
df <- read.table("pred_svm_sep.txt", header = T)
plot_ly(data = df, x = ~obs, y = ~pred, text=~line) %>% add_markers(alpha = 0.5)

x_mid <- mean(c(max(df$obs, na.rm = T),
                min(df$obs, na.rm = T)))
y_mid <- mean(c(max(df$pred, na.rm = T),
                min(df$pred, na.rm = T)))

df %>% mutate(quadrant = case_when(obs > x_mid & pred > y_mid ~ "Q1",
                                   obs <= x_mid & pred > y_mid ~ "Q2",
                                   obs <= x_mid & pred <= y_mid ~ "Q3",
                                   TRUE ~ "Q4")) %>% 
  ggplot(aes(x = obs, y = pred, color = quadrant)) + geom_vline(xintercept = x_mid) + # plot vertical line
  geom_hline(yintercept = y_mid) + # plot horizontal line
  geom_point() + theme_classic(base_size = 16) + theme(legend.position = "none")


library(dplyr)
library(ggplot2)
x_mid <- mean(c(max(iris$Petal.Length, na.rm = TRUE), 
                min(iris$Petal.Length, na.rm = TRUE)))

y_mid <- mean(c(max(iris$Petal.Width, na.rm = TRUE), 
                min(iris$Petal.Width, na.rm = TRUE)))



iris %>% 
  mutate(quadrant = case_when(Petal.Length > x_mid & Petal.Width > y_mid   ~ "Q1",
                              Petal.Length <= x_mid & Petal.Width > y_mid  ~ "Q2",
                              Petal.Length <= x_mid & Petal.Width <= y_mid ~ "Q3",
                              TRUE                                         ~ "Q4")) %>% 
  ggplot(aes(x = Petal.Length, y = Petal.Width, color = quadrant)) +
  geom_vline(xintercept = x_mid) + # plot vertical line
  geom_hline(yintercept = y_mid) + # plot horizontal line
  geom_point()




plot_ly(data = df, x = ~obs, y = ~pred, mode = "markers+text", text= ~line,
        textposition = "top", hoverinfo = "x+y") %>% add_markers(alpha = 0.5, ids = ~line) %>% highlight(on = "plotly_hover", off = "plotly_doubleclick")




