#####define_function##############################################################
l=length; co=colors();
co_hsa=co[55];  co_ptr=co[81];  co_mml=co[30];
#########################
age.test.polynomials.1.f = function(Mat, Ages, verbose=F, pvals=FALSE, PCUTOFF=0.05) {
#tries different models using one variable up to 3rd degree. Returns the index of the best model in Models[[3]]
   parI = list(Ages, Ages^2, Ages^3) 
   Genes = nrow(Mat); PI = length(parI); LI = 2^ (PI); Models = Models[[3]]
   LmListInt1 = function (y, par, m) 
   { lm(y ~ I(par[[1]] * m[1]) + I(par[[2]] * m[2]) + I(par[[3]] * m[3])) }
   if (pvals==FALSE) {
      Result = apply(Mat, 1, function(y) {
         lr = lapply(1:8, function(j) {LmListInt1(y, parI, Models[j,])})
         mod = sapply(2:8, function(k) {anova(lr[[1]], lr[[k]])$Pr[2]})
         mod[is.na(mod)] = 1
         ifelse(min(mod) <= PCUTOFF, which.min(mod)+1, 1)
      })
      names(Result) = rownames(Mat)
   } else {
      Result = t(apply(Mat, 1, function(y) {
         lr = lapply(1:8, function(j) {LmListInt1(y, parI, Models[j,])})
         mod = sapply(2:8, function(k) {anova(lr[[1]], lr[[k]])$Pr[2]})
         mod[is.na(mod)] = 1
         c(ifelse(min(mod) <= PCUTOFF, which.min(mod)+1, 1), min(mod))
      }))
      rownames(Result) = rownames(Mat)
      colnames(Result) = c("Model", "P-value")
   }
   return(Result)
}
#########################
Models = list()
Models[[6]] = matrix(0, 64, 6)
{ j = 1     
{
for (i1 in 0:1) {
   for (i2 in 0:1) {
      for (i3 in 0:1) {
         for (i4 in 0:1) {
            for (i5 in 0:1) {
               for (i6 in 0:1) {
               Models[[6]][j,] = c(i1,i2,i3,i4,i5,i6)
               j = j + 1 } } } } } } } }

Models[[5]] = matrix(0, 32, 5)
{ j = 1     
{
for (i1 in 0:1) {
   for (i2 in 0:1) {
      for (i3 in 0:1) {
         for (i4 in 0:1) {
            for (i5 in 0:1) {
               Models[[5]][j,] = c(i1,i2,i3,i4, i5)
               j = j + 1 } } } } } } }
               
Models[[4]] = matrix(0, 16, 4)
{ j = 1     
{
for (i1 in 0:1) {
   for (i2 in 0:1) {
      for (i3 in 0:1) {
         for (i4 in 0:1) {
            Models[[4]][j,] = c(i1,i2,i3,i4)
               j = j + 1 } } } } } }

Models[[3]] = matrix(0, 8, 3)
{ j = 1     
{
for (i1 in 0:1) {
   for (i2 in 0:1) {
      for (i3 in 0:1) {
         Models[[3]][j,] = c(i1,i2,i3)
               j = j + 1 } } } } }
                              
Models[[2]] = matrix(0, 4, 2)
{ j = 1     
{
for (i1 in 0:1) {
   for (i2 in 0:1) {
      Models[[2]][j,] = c(i1,i2)
               j = j + 1 } } } }

Models[[1]] = matrix(0, 2, 1)
{ j = 1     
{
for (i1 in 0:1) {
   Models[[1]][j,] = c(i1)
               j = j + 1 } } }              
##############################
#function to calculate proportion of variance explained by age
prop.var.age.1.f = 
function(Mat, Ages) {
  x0 = apply(Mat, 1, function(y) { tryCatch(var(y)*(length(y)-1), error = function(e) NA) })  #null = total variance = age effect + biological noise + technical noise
  x1 = apply(Mat, 1, function(y) { 
    x = Ages
    m2 = lm(y ~ I(x) * 1 + I(x^2) * 1 + I(x^3) * 1)
    tryCatch(sum(m2$resi^2), error = function(e) NA) })    #only age, degree=3. biological noise + technical noise.
  100*(1 - x1/x0)   #age effect
}
##############################
normalize.rows.f = function (Mat) {
    if (length(nrow(Mat)) != 0) {
        t(apply(Mat, 1, function(y) (y-mean(y))/sd(y))) 
    } else {
        c(apply(t(Mat), 1, function(y) (y-mean(y))/sd(y))) 
    }
}

#FUNCTIONS FOR TESTING DIFFERENTIAL EXPRESSION 
species.test.multiregression.1.f = function(
#F-test for species difference, where the regression curve is chosen for refspecies. 
   E,           #expression matrix
   x,           #ages
   RefSpecies,  #indicator of species 1
   Species2,    #indicator of species 2
   DEGREE = 3,            #degree of the polynomial regression model to use. obsolete. 
   EACH = TRUE,           #if EACH is TRUE a list of submodels are tested for species difference, where the second species can have its own parameters. If EACH is FALSE, only one model is tested, where the second species has its own parameters for each element of the model.
   Ready.Bm = c(),         #vector of previously calculated best expression-age regression models. if left empty, calculates the models de novo. 
   SIG.MODELS = TRUE       #should we only use significant expression-age models, or just the best model irrespective of significance
) {
   mff = age.test.polynomials.1.f; mo1 = multireg.models.1.f; 
   	if (EACH) mo2 = ancova.models.02.f else mo2 = ancova.models.01.f

	 if (l(Ready.Bm) != 0) {					#find best model or get it
		 BmMat = Ready.Bm
	 } else { BmMat = mff(E[,RefSpecies], x[RefSpecies], verbose=F, PCUTOFF=ifelse(SIG.MODELS, 0.05, 1)) }									

   Sp = factor(c(rep(1, l(RefSpecies)), rep(0, l(Species2))))
   XYZ = c(RefSpecies, Species2)
   Res = t(sapply(1:nrow(E), function(i) {
      y = E[i,]
      Bm = BmMat[i]                                             
      null = mo1(Bm, x[XYZ], y[XYZ])
      altlist = mo2(Bm, x[XYZ], y[XYZ], Sp)
      if (EACH & Bm != 1) pval = sapply(altlist, function(alt) anova(alt, null)$Pr[2]) else pval = anova(altlist, null)$Pr[2]
  #########################################    
      if(min(pval) == "NaN") c(1, Bm, 1) else c(min(pval), Bm, which.min(pval))
  #########################################
      #c(min(pval), Bm, which.min(pval))
   }))
   rownames(Res) = rownames(E)
   colnames(Res) = c("P-value", "Multireg_Model", "Ancova_Model")
   return(Res)
}

###############################
species.test.multiregression.1.f2 = function(
#F-test for species difference, where the regression curve is chosen for refspecies. 
   E,           #expression matrix
   x,           #ages
   RefSpecies,  #indicator of species 1
   Species2,    #indicator of species 2
   DEGREE = 3,            #degree of the polynomial regression model to use. obsolete. 
   EACH = TRUE,           #if EACH is TRUE a list of submodels are tested for species difference, where the second species can have its own parameters. If EACH is FALSE, only one model is tested, where the second species has its own parameters for each element of the model.
   Ready.Bm = c(),         #vector of previously calculated best expression-age regression models. if left empty, calculates the models de novo. 
   SIG.MODELS = TRUE       #should we only use significant expression-age models, or just the best model irrespective of significance
) {
   mff = age.test.polynomials.1.f; mo1 = multireg.models.1.f; 
   	if (EACH) mo2 = ancova.models.02.f else mo2 = ancova.models.01.f

	 if (l(Ready.Bm) != 0) {					#find best model or get it
		 BmMat = Ready.Bm
	 } else { BmMat = mff(E[,RefSpecies], x[RefSpecies], verbose=F, PCUTOFF=ifelse(SIG.MODELS, 0.05, 1)) }									

   Sp = factor(c(rep(1, l(RefSpecies)), rep(0, l(Species2))))
   XYZ = c(RefSpecies, Species2)
   Res = t(sapply(1:nrow(E), function(i) {
      y = E[i,]
      Bm = BmMat[i]                                             
      null = mo1(Bm, x[XYZ], y[XYZ])
      altlist = mo2(Bm, x[XYZ], y[XYZ], Sp)
      if (EACH & Bm != 1) pval = sapply(altlist, function(alt) anova(alt, null)$Pr[2]) else pval = anova(altlist, null)$Pr[2]
  #########################################    
      #if(min(pval) == "NaN") c(1, Bm, 1) else c(min(pval), Bm, which.min(pval))
  #########################################
      test = c(min(pval), Bm, which.min(pval));
      if (length(test) == 2) c(1,Bm,1) else (c(min(pval), Bm, which.min(pval)))
   }))
   rownames(Res) = rownames(E)
   colnames(Res) = c("P-value", "Multireg_Model", "Ancova_Model")
   return(Res)
}
#####################################

species.test.multiregression.1.f3 = function(
#F-test for species difference, where the regression curve is chosen for refspecies. 
   E,           #expression matrix
   x,           #ages
   RefSpecies,  #indicator of species 1
   Species2,    #indicator of species 2
   DEGREE = 3,            #degree of the polynomial regression model to use. obsolete. 
   EACH = TRUE,           #if EACH is TRUE a list of submodels are tested for species difference, where the second species can have its own parameters. If EACH is FALSE, only one model is tested, where the second species has its own parameters for each element of the model.
   Ready.Bm = c(),         #vector of previously calculated best expression-age regression models. if left empty, calculates the models de novo. 
   SIG.MODELS = TRUE       #should we only use significant expression-age models, or just the best model irrespective of significance
) {
   mff = age.test.polynomials.1.f; mo1 = multireg.models.1.f; 
   	if (EACH) mo2 = ancova.models.02.f else mo2 = ancova.models.01.f

	 if (l(Ready.Bm) != 0) {					#find best model or get it
		 BmMat = Ready.Bm
	 } else { BmMat = mff(E[,RefSpecies], x[RefSpecies], verbose=F, PCUTOFF=ifelse(SIG.MODELS, 0.05, 1)) }									

   Sp = factor(c(rep(1, l(RefSpecies)), rep(0, l(Species2))))
   XYZ = c(RefSpecies, Species2)
   Res = t(sapply(1:nrow(E), function(i) {
      y = E[i,]
      Bm = BmMat[i]                                             
      null = mo1(Bm, x[XYZ], y[XYZ])
      altlist = mo2(Bm, x[XYZ], y[XYZ], Sp)
      if (EACH & Bm != 1) pval = sapply(altlist, function(alt) anova(alt, null)$Pr[2]) else pval = anova(altlist, null)$Pr[2]
  #########################################    
      #if(min(pval) == "NaN") c(1, Bm, 1) else c(min(pval), Bm, which.min(pval))
  #########################################
      #test = c(min(pval), Bm, which.min(pval));
      #if (length(test) == 2) c(1,Bm,1) else (c(min(pval), Bm, which.min(pval)))
  #########################################
        for (n in 1:length(pval)) {
        if (is.na(pval[n]))  {pval[n] = 1;}
        if (pval[n] == "NaN")  {pval[n] = 1;}
        }
        c(min(pval), Bm, which.min(pval));
   }))
   rownames(Res) = rownames(E)
   colnames(Res) = c("P-value", "Multireg_Model", "Ancova_Model")
   return(Res)
}
#####################################

multireg.models.1.f = function(BestModels, x, y)  {
   if (BestModels == 1) m2 = lm(y ~ rep(1, l(x)) * 1) 
   if (BestModels == 2) m2 = lm(y ~ I(x^3) * 1)
   if (BestModels == 3) m2 = lm(y ~ I(x^2) * 1)
   if (BestModels == 4) m2 = lm(y ~ I(x^2) * 1 + I(x^3) * 1)
   if (BestModels == 5) m2 = lm(y ~ I(x) * 1)
   if (BestModels == 6) m2 = lm(y ~ I(x) * 1 + I(x^3) * 1)
   if (BestModels == 7) m2 = lm(y ~ I(x) * 1 + I(x^2) * 1)
   if (BestModels == 8) m2 = lm(y ~ I(x) * 1 + I(x^2) * 1 + I(x^3) * 1)
   m2
}

ancova.models.01.f = function(BestModels, x, y, Sp)  {
#each species has specific slope and intercept
   if (BestModels == 1) m2 = lm(y ~ rep(1, l(x)) * Sp) 
   if (BestModels == 2) m2 = lm(y ~ I(x^3) * Sp)
   if (BestModels == 3) m2 = lm(y ~ I(x^2) * Sp)
   if (BestModels == 4) m2 = lm(y ~ I(x^2) * Sp + I(x^3) * Sp)
   if (BestModels == 5) m2 = lm(y ~ I(x) * Sp)
   if (BestModels == 6) m2 = lm(y ~ I(x) * Sp + I(x^3) * Sp)
   if (BestModels == 7) m2 = lm(y ~ I(x) * Sp + I(x^2) * Sp)
   if (BestModels == 8) m2 = lm(y ~ I(x) * Sp + I(x^2) * Sp + I(x^3) * Sp)
   m2
}

ancova.models.02.f = function(BestModels, x, y, Sp)  {
#a list of models where each species *may have* a specific slope and intercept
   Sp2 = rep(1, l(Sp))
   if (BestModels == 1) 
      m2 = lm(y ~ rep(1, l(x)) * Sp) 
   if (BestModels == 2) 
      m2 = lapply(1:2, function(h) {modlist = list(Sp2,Sp)[(Models[[1]]+1)[h,]]; lm(y ~ I(x^3) * modlist[[1]] + Sp)})
   if (BestModels == 3) 
      m2 = lapply(1:2, function(h) {modlist = list(Sp2,Sp)[(Models[[1]]+1)[h,]]; lm(y ~ I(x^2) * modlist[[1]] + Sp)})
   if (BestModels == 4) 
      m2 = lapply(1:4, function(h) {modlist = list(Sp2,Sp)[(Models[[2]]+1)[h,]]; lm(y ~ I(x^2) * modlist[[1]] + I(x^3) * modlist[[2]] + Sp)})
   if (BestModels == 5) 
      m2 = lapply(1:2, function(h) {modlist = list(Sp2,Sp)[(Models[[1]]+1)[h,]]; lm(y ~ I(x) * modlist[[1]] + Sp)})
   if (BestModels == 6) 
      m2 = lapply(1:4, function(h) {modlist = list(Sp2,Sp)[(Models[[2]]+1)[h,]]; lm(y ~ I(x) * modlist[[1]] + I(x^3) * modlist[[2]] + Sp)})
   if (BestModels == 7) 
      m2 = lapply(1:4, function(h) {modlist = list(Sp2,Sp)[(Models[[2]]+1)[h,]]; lm(y ~ I(x) * modlist[[1]] + I(x^2) * modlist[[2]] + Sp)})
   if (BestModels == 8) 
      m2 = lapply(1:8, function(h) {modlist = list(Sp2,Sp)[(Models[[3]]+1)[h,]]; lm(y ~ I(x) * modlist[[1]] + I(x^2) * modlist[[2]] + I(x^3) * modlist[[3]] + Sp)})
   m2
}

k2try.1.f = function(
  Mat,  #normalized expression matrix
  k2try,  #vector of integers indicating which K's to try (eg. c(4,6,8))
  REPEAT=1000   #number of repetitions
) {
    km_l1 = lapply(k2try, function(k) { print(k)
            tab = sapply(1:REPEAT, function(i) {sort(table(kmeans(Mat, k)$cluster))} );
            tabx = rev(sort(table(apply(tab, 2, function(x) paste(x, collapse="_"))))); tabx }); names(km_l1) = k2try
    Tab = sapply(km_l1, function(x) as.numeric(x[1:3])); print(Tab)
    return(km_l1)
}

postk2try.1.f = function(
  Mat,  #normalized expression matrix 
  KM,   #integer indicating which K to use
  List  #output of k2try.1.f
) {
    km1 = { 
            i = 0; while (i != 1) {
                    cat(i); km1 = kmeans(Mat, KM);
                    if (paste(sort(table(km1$cluster)), collapse="_") == names(List[[as.character(KM)]][1])) { KmX = km1; i = 1 } }
      Kmc = KmX$cluster; Kmc2 = Kmc; 
      for (J in 1:KM) { Kmc2[Kmc %in% names(rev(sort(table(Kmc))))[J]] = J }
      KmX$cluster = Kmc2; KmX } 
  return(km1)
}
