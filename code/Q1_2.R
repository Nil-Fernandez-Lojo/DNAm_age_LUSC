rm(list = ls())

library(survival)
library(car)
library(rstatix)
alpha = 0.05
sidak <- function(alpha, N){
    return(1-(1-alpha)^(1/N))
}
#Load DNAm data (variable name dm.age.subset)
load('../data/LUSC-1.Rda')

#Load Clinical data
df <- as.data.frame(t(read.table(file = '../data/Clinical/LUSC.clin.merged.picked.txt', sep = '\t', header = TRUE,row.names=1)))
index <- c("years_to_birth", "days_to_death", "days_to_last_followup", "date_of_initial_pathologic_diagnosis", "days_to_last_known_alive", "karnofsky_performance_score", "number_pack_years_smoked", "year_of_tobacco_smoking_onset")
df[index] <- sapply(df[index], function(x) as.numeric(as.character(x)))

#process smoking variables
df$number_years_smoked <- df$date_of_initial_pathologic_diagnosis - df$year_of_tobacco_smoking_onset
df$total_number_packs_smoked <- df$number_years_smoked * df$number_pack_years_smoked

#Use NA notation instead of x notation
df$pathology_M_stage[df$pathology_M_stage=="mx"]<-NA
df$pathology_N_stage[df$pathology_N_stage=="nx"]<-NA
df$residual_tumor[df$residual_tumor=="rx"]<-NA
df <- droplevels(df)

#Merge some cancer stages (e.g. ii,iia,iib into single stage 2) and convert them into numerical value
#pathologic stage
df$pathologic_stage_merged <- df$pathologic_stage
levels(df$pathologic_stage_merged) <- c("1","1","1","2","2","2", "3", "3", "3","4")
df$pathologic_stage_numeric <-  as.numeric(df$pathologic_stage)
df$pathologic_stage_merged_numeric <-  as.numeric(df$pathologic_stage_merged)

#T stage
df$pathology_T_stage_merged <- df$pathology_T_stage
levels(df$pathology_T_stage_merged) <- c("1", "1","1", "2", "2","2","3","4")
df$pathology_T_stage_merged_numeric <-  as.numeric(df$pathology_T_stage_merged)

#M stage
df$pathology_M_stage_merged <- df$pathology_M_stage
levels(df$pathology_M_stage_merged) <- c("0", "1","1", "1")
df$pathology_M_stage_merged_numeric <-  as.numeric(df$pathology_M_stage_merged)

#N stage
df$pathology_N_stage_numeric <-  as.numeric(df$pathology_N_stage)

#add TSS variable and Participant to Clinical data
TCGA_names_df_split <- data.frame(do.call(rbind,strsplit(rownames(df), ".", fixed=TRUE))[,3:2])
colnames(TCGA_names_df_split) <- c("Participant","TSS")
df <- cbind(TCGA_names_df_split,df)

#generate datframe with DNAm informtation and split it according to sample type
DNAm_df <- data.frame(do.call(rbind,strsplit(names(dm.age.subset), ".", fixed=TRUE))[,3:4])
colnames(DNAm_df) <- c("Participant","Sample")
levels(DNAm_df$Sample) <- c("tumor", "tumor","healthy", "healthy")
DNAm_df$DNAm_age <- dm.age.subset
DNAm_df_split<-split(DNAm_df, DNAm_df$Sample)
DNAm_healthy <-subset(DNAm_df_split$healthy, select = -Sample)
names(DNAm_healthy)[names(DNAm_healthy) == 'DNAm_age'] <- 'DNAm_age_healthy'
DNAm_tumor <-subset(DNAm_df_split$tumor, select = -Sample)
names(DNAm_tumor)[names(DNAm_tumor) == 'DNAm_age'] <- 'DNAm_age_tumor'

# using any(table(DNAm_healthy$Participant)>1) and any(table(DNAm_tumor$Participant)>1) we can check that unique values per participant

#merge DNAm_df with df
df <- merge(df, DNAm_healthy, by = "Participant",all = TRUE)
df <- merge(df, DNAm_tumor, by = "Participant",all = TRUE)
linearMod <- lm(DNAm_age_healthy ~ years_to_birth, data=df)
#distPred <- predict(lmMod, df$years_to_birth[!is.na()])  # predict distance

newdf = data.frame(years_to_birth = df$years_to_birth[!is.na(df$DNAm_age_tumor)])
df$DNAm_age_acceleration <- df$DNAm_age_tumor - predict(linearMod, df)

Q1_tests = c("m age healthy - age", "m age tumor - age", "m age tumor - m age healthy", "m acc - m age healthy")
Q1_corr = data.frame(p = numeric(length(Q1_tests)), rho = numeric(length(Q1_tests)), significant =  numeric(length(Q1_tests)))
rownames(Q1_corr) = Q1_tests
a  <- cor.test(x=df$years_to_birth, y= df$DNAm_age_healthy, method = 'spearman')
Q1_corr[1,"p"] <- a[[3]]
Q1_corr[1,"rho"] <- a[[4]]
Q1_1_pearson <- cor.test(x=df$years_to_birth, y= df$DNAm_age_healthy, method = 'pearson') 
a <- cor.test(x=df$years_to_birth, y= df$DNAm_age_tumor, method = 'spearman')
Q1_corr[2,"p"] <- a[[3]]
Q1_corr[2,"rho"] <- a[[4]]
a <- cor.test(x=df$DNAm_age_healthy, y= df$DNAm_age_tumor, method = 'spearman')
Q1_corr[3,"p"] <- a[[3]]
Q1_corr[3,"rho"] <- a[[4]]
a <- cor.test(x=df$years_to_birth, y= df$DNAm_age_acceleration, method = 'spearman')
Q1_corr[4,"p"] <- a[[3]]
Q1_corr[4,"rho"] <- a[[4]]

Q1_corr$significant <- Q1_corr$p <= sidak(alpha,length(Q1_tests))

#ANOVA for patients for which we have the DNAm age for both healthy and cancer tissues 
DNAm_age <- c(df$DNAm_age_healthy[!is.na(df$DNAm_age_healthy + df$DNAm_age_tumor)], df$DNAm_age_tumor[!is.na(df$DNAm_age_healthy + df$DNAm_age_tumor)])
Chron_age <- c(df$years_to_birth[!is.na(df$DNAm_age_healthy + df$DNAm_age_tumor)], df$years_to_birth[!is.na(df$DNAm_age_healthy + df$DNAm_age_tumor)])
group <- as.factor(c(rep("Healthy",sum(!is.na(df$DNAm_age_healthy + df$DNAm_age_tumor))),rep("Tumour",sum(!is.na(df$DNAm_age_healthy + df$DNAm_age_tumor)))))

lin_model_2 <- lm(DNAm_age ~ Chron_age+group)
summary(lin_model_2)

# Q2 correlations
var_pearson_corr <- c("karnofsky_performance_score",
                      "number_years_smoked", 
                      "total_number_packs_smoked", 
                      "number_pack_years_smoked")

var_cat = c("pathologic_stage_merged",
            "pathology_T_stage_merged",
            "pathology_M_stage_merged",
            "pathology_N_stage",
            "gender",
            "radiation_therapy",
            "histological_type",
            "residual_tumor",
            "ethnicity",
            "race"
)

#lin_model_3 <- lm(DNAm_age_tumor ~ karnofsky_performance_score
#                  +number_years_smoked
#                  +total_number_packs_smoked
#                  +number_pack_years_smoked
#                  +pathology_T_stage_merged
#                  +pathology_M_stage_merged
#                  +pathology_N_stage
#                  +gender
#                  +radiation_therapy
#                  +histological_type
#                  +residual_tumor
#                  +ethnicity
#                  +race
#                  data = df)

corr_pearson <- data.frame(rho_age = numeric(length(var_pearson_corr)),
                           p.value_age= numeric(length(var_pearson_corr)),
                           rho_acc = numeric(length(var_pearson_corr)),
                           p.value_acc= numeric(length(var_pearson_corr)))
rownames(corr_pearson) <- var_pearson_corr
for (var_corr in var_pearson_corr) {
    a <- cor.test(x=df[[var_corr]], y= df$DNAm_age_tumor, method = 'spearman')
    corr_pearson[var_corr,"rho_age"] <- a$estimate
    corr_pearson[var_corr, "p.value_age"] <- a$p.value
    a <- cor.test(x=df[[var_corr]], y= df$DNAm_age_acceleration, method = 'spearman')
    corr_pearson[var_corr,"rho_acc"] <- a$estimate
    corr_pearson[var_corr, "p.value_acc"] <- a$p.value
}
corr_pearson$significance_age <- corr_pearson$p.value_age <= sidak(alpha, length(var_pearson_corr)+length(var_cat))
corr_pearson$significance_acc <- corr_pearson$p.value_acc <= sidak(alpha, length(var_pearson_corr)+length(var_cat))

df_categorical = df[var_cat]
corr_cat <- data.frame(p_levene_age = rep(NA,length(var_cat)),
                       levene_age_passed = NA,
                       p_ANOVA_age = NA,
                       ANOVA_age_significant =NA,
                       p_shapiro_age = NA,
                       shapiro_age_passed = NA,
                       p_levene_acc = NA,
                       p_Welch_ANOVA_age = NA,
                       Welch_ANOVA_age_significant = NA,
                       p_kruskal_age =NA,
                       kruskal_age_significant =NA,
                       levene_acc_passed = NA,
                       p_shapiro_acc = NA,
                       shapiro_acc_passed = NA,
                       p_ANOVA_acc = NA,
                       ANOVA_acc_significant = NA,
                       p_kruskal_acc = NA,
                       kruskal_acc_significant = NA)
rownames(corr_cat) <- var_cat

corr_cat$p_levene_age <- apply(df_categorical,2,function(x) {leveneTest(DNAm_age_tumor ~ x, data = df)[1,3]})
corr_cat$p_levene_acc <- apply(df_categorical,2,function(x) {leveneTest(DNAm_age_acceleration ~ x, data = df)[1,3]})
corr_cat$levene_age_passed <- corr_cat$p_levene_age > alpha
corr_cat$levene_acc_passed <- corr_cat$p_levene_acc > alpha

res_aov_age <- list()
res_aov_acc <- list()
Normality_Welch <- list()
for (i in 1:length(var_cat)) {
    if (corr_cat$levene_age_passed[i]) {
        res_aov_age[[i]] <- aov(DNAm_age_tumor ~ df[[var_cat[i]]], data = df)
        a = summary(res_aov_age[[i]])
        corr_cat[i,"p_ANOVA_age"] <-a[[1]][1,'Pr(>F)']
        aov_residuals <- residuals(object = res_aov_age[[i]])
        corr_cat[i,"p_shapiro_age"] <-shapiro.test(x = aov_residuals )[2]
    }   
    else {
        res_aov_age[[i]] <- oneway.test(DNAm_age_tumor ~ df[[var_cat[i]]], data = df)
        corr_cat[i,"p_Welch_ANOVA_age"] <- res_aov_age[[i]][3]
        a <- list()
        for (j in 1:length(levels(df[[var_cat[i]]][!is.na(df$DNAm_age_tumor)]))) {
            a[j] <- shapiro.test(df$DNAm_age_tumor[!is.na(df$DNAm_age_tumor) & (df[[var_cat[i]]] == levels(df[[var_cat[i]]])[j])])[2]
        }
        Normality_Welch <- append(Normality_Welch,a)
    }
    res_aov_acc[[i]] <- aov(DNAm_age_acceleration ~ df[[var_cat[i]]], data = df)
    a = summary(res_aov_acc[[i]])
    corr_cat[i,"p_ANOVA_acc"] <-a[[1]][1,'Pr(>F)']
    aov_residuals <- residuals(object = res_aov_acc[[i]])
    corr_cat[i,"p_shapiro_acc"] <-shapiro.test(x = aov_residuals )[2]
}
corr_cat$Welch_ANOVA_age_significant <- corr_cat$p_Welch_ANOVA_age < sidak(alpha, length(var_cat)+length(var_pearson_corr))
corr_cat$ANOVA_age_significant <- corr_cat$p_ANOVA_age < sidak(alpha, length(var_cat)+length(var_pearson_corr))
corr_cat$ANOVA_acc_significant <- corr_cat$p_ANOVA_acc < sidak(alpha, length(var_cat)+length(var_pearson_corr))
corr_cat$shapiro_age_passed <- corr_cat$p_shapiro_age > alpha
corr_cat$shapiro_acc_passed <- corr_cat$p_shapiro_acc >alpha

var_kruskal_age = na.omit(var_cat[!corr_cat$shapiro_age_passed])
res_kruskal_age <- list()
for (i in 1:length(var_kruskal_age)) {
    res_kruskal_age[[i]] <- kruskal.test(DNAm_age_tumor ~ df[[var_kruskal_age[i]]], data = df)
    corr_cat[var_kruskal_age[i],"p_kruskal_age"] <- res_kruskal_age[[i]]$p.value
}

var_kruskal_acc = var_cat[!corr_cat$shapiro_acc_passed]
res_kruskal_acc <- list()
for (i in 1:length(var_kruskal_acc)) {
    res_kruskal_acc[[i]] <- kruskal.test(DNAm_age_acceleration ~ df[[var_kruskal_acc[i]]], data = df)
    corr_cat[var_kruskal_acc[i],"p_kruskal_acc"] <- res_kruskal_acc[[i]]$p.value
}

corr_cat$kruskal_age_significant <- corr_cat$p_kruskal_age < sidak(alpha, length(var_cat)+length(var_pearson_corr))
corr_cat$kruskal_acc_significant <- corr_cat$p_kruskal_acc < sidak(alpha, length(var_cat)+length(var_pearson_corr))

# Only race statistically significant (but probably due to 1 outlier sample 98, participant 3407)
df_remove_outlier <- df[-98,]
aov_res_race_rem_outlier <- aov(DNAm_age_acceleration ~ race, data = df_remove_outlier)


#Survival analysis
df$DNAm_age_tumor_high = as.factor(df$DNAm_age_tumor > mean(df$DNAm_age_tumor, na.rm=TRUE))
df$DNAm_acc_tumor_high = as.factor(df$DNAm_age_acceleration > mean(df$DNAm_age_acceleration, na.rm=TRUE))

ev = as.numeric(df$vital_status)
time_event <- df$days_to_death
time_event[is.na(time_event)] <- df$days_to_last_followup[is.na(time_event)]
su = Surv(time_event, ev)
#km_age = survfit(su~DNAm_age_tumor_high, data=df)
km_acc = survfit(su~DNAm_acc_tumor_high, data=df)

res_cox_age <- coxph(su ~ DNAm_age_tumor, data =  df)
res_cox_acc <- coxph(su ~ DNAm_age_acceleration, data =  df)

#Graph generation:
# Figure 1
png(file="../report/img/fig1.png",
    units="in", 
    width=10, 
    height=5,
    res=300)
plot(df$years_to_birth, df$DNAm_age_healthy,
     xlab = "Chronological age (years)",
     ylab = "DNAm age",
)
abline(linearMod, col = "red")
abline(0,1, col = "blue")
grid()
title(main= paste("DNAm age healthy tissue vs Chronological age \n","Spearman r =", as.character(signif(Q1_corr[1,"rho"],3)), "p = ", as.character(signif(Q1_corr[1,"p"],3))))
legend("topleft", 
       legend=c("Data", "Linear regression", "y=x"),
       col=c("black", "red", "blue"),
       pch = c(1,NA,NA),
       lty =c(NA,1,1))
dev.off()

# Figure 2
png(file="../report/img/fig2.png",
    units="in", 
    width=10, 
    height=5,
    res=300)
plot(df$years_to_birth, df$DNAm_age_tumor,
     xlab = "Chronological age (years)",
     ylab = "DNAm age",
)
abline(linearMod, col = "red")
grid()
title(main=paste("DNAm age tumour tissue vs Chronological age \n", "Spearman r =", as.character(signif(Q1_corr[2,"rho"],3)), "p = ", as.character(signif(Q1_corr[2,"p"],3))))
legend("topleft", 
       legend=c("Data", "Linear regression prediction"),
       col=c("black", "red"),
       pch = c(1,NA),
       lty =c(NA,1))
dev.off()

# Figure 3
png(file="../report/img/fig3.png",
    units="in", 
    width=10, 
    height=5,
    res=300)
plot(df$DNAm_age_healthy, df$DNAm_age_tumor,
     xlab = "Healthy tissue DNAm age",
     ylab = "Tumour tissue DNAm age",
)
abline(0,1, col = "blue")
grid()
title(main=paste("Tumour vs healthy tissue DNAm age \n","Spearman r =", as.character(signif(Q1_corr[3,"rho"],3)), "p = ", as.character(signif(Q1_corr[3,"p"],3))))
legend("topleft", 
       legend=c("Data", "y=x"),
       col=c("black", "blue"),
       pch = c(1,NA),
       lty =c(NA,1))
dev.off()

# Figure 4
png(file="../report/img/fig4.png",
    units="in", 
    width=10, 
    height=5,
    res=300)
plot(df$years_to_birth, df$DNAm_age_acceleration,
     xlab = "DNAm acceleration",
     ylab = "Chronological age",
)
abline(0,1, col = "blue")
grid()
title(main=paste("Tumour DNAm acceleration vs chronological age \n","Spearman r =", as.character(signif(Q1_corr[4,"rho"],3)), "p = ", as.character(signif(Q1_corr[4,"p"],3))))

dev.off()

# Figure 5
png(file="../report/img/fig5.png",
    units="in", 
    width=10, 
    height=5,
    res=300)
plot(km_acc,
     col = c("red", "blue"),
     xlab = "Time (days)",
     ylab = "Survival probability")
legend("topright", 
       legend=c("High DNAm acceleration", "Low DNAm acceleration"),
       col=c("red", "blue"),
       lty =c(1,1))
grid()
title(main = "Kaplan–Meier plot of cancer population survival")
dev.off()
