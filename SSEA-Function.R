# Funktionen, um SSEA durchzuführen

### Ein paar Zufallszahlen
p.rand <- matrix(round(runif(2600, 0, 1), 8), ncol=26, byrow=T)
p.obs <- matrix(round(runif(26, 0, 1), 8), ncol=26, byrow=T)
colnames(p.rand) <- letters[1:26]
colnames(p.obs) <- letters[1:26]


set1 <- letters[c(1:4, 6)]
set2 <- letters[c(5,7:13)]
set3 <- letters[14:25]
set4 <- letters[26]

liste.obs <- vector(mode="list", length=4)
liste.obs[[1]] <- p.obs[,colnames(p.obs) %in% set1]
liste.obs[[2]] <- p.obs[,colnames(p.obs) %in% set2]
liste.obs[[3]] <- p.obs[,colnames(p.obs) %in% set3]
liste.obs[[4]] <- p.obs[,colnames(p.obs) %in% set4]
names(liste.obs) <- c("set1", "set2", "set3", "set4")


liste.rand <- vector(mode="list", length=4)
liste.rand[[1]] <- p.rand[,colnames(p.rand) %in% set1]
liste.rand[[2]] <- p.rand[,colnames(p.rand) %in% set2]
liste.rand[[3]] <- p.rand[,colnames(p.rand) %in% set3]
liste.rand[[4]] <- p.rand[,colnames(p.rand) %in% set4]
names(liste.rand) <- c("set1", "set2", "set3", "set4")




W_k_perm <- function(perm.data, k_first){
  #~ auslagerbar
  is.wholenumber <- function(k_first, tol = .Machine$double.eps^0.5){ ### Abfangfunktion, ob k_first
    abs(k_first - round(k_first)) < tol}                                # ganzzahlig ist
  #~
  if(is.list(perm.data)!=TRUE){                                           ### Fehlerabfangung; Input muss eine 
                                                                        # Liste sein; jedes Element der Liste 
                                                                        # ist ein Dataframe/matrix mit SNPs
    stop("Stop! Input has to be a list. Each element hast to be a SNP-data-frame for a specific gen")
  }
  if(is.wholenumber(k_first)!=TRUE | k_first<=0){
    stop("Stop! k_first has to be a whole positiv number")
  }
  
  
  else{
    lapply(c(seq_along(perm.data)), FUN=function(var){   ### äußere lapply für Genauswahl der in der Liste
      
      if(length(dim(perm.data[[var]]))>1){      ### Abfrage, ob Gen nicht nur aus einem SNP besteht; dann
                                                  # muss apply nicht angewendet werden
        apply(perm.data[[var]], MARGIN=1,             ### inneres apply für zeilenweisen Zugriff auf die 
                                                        # Datensätze/matrizen mit den Variablen
              
              FUN=function(x){
                sorted.p <- sort(x)[1:k_first]         ### Sortierung der p.-Werte und anschließend
                sorted.p <- sorted.p[!is.na(sorted.p)]   # nur die behalten, die auch wirkluch da sind (NAs raus)
                p.product <- prod(sorted.p)
              } #((())) Ende der Funktion vom inneren apply
              
        ) #((())) Ende inneres apply
      } #((())) Ende inneres if
      
      else{   ### wird ausgewählt, wenn nur ein SNP in einem Gen vorliegt
        p.product <- perm.data[[var]] ### Wenn es nur einen SNP gibt, muss nichts sortiert/multipliziert werden
      } #((())) Ende inneres else
      
      
    } #((())) Ende function von äußerem lapply
    
    ) #((())) Ende äußeres lapply
  
  } #((())) Ende von else
  
} #((())) Ende der Funktion



# für permutierte daten; nur zum späteren Überprüfen
per_k1 <- W_k_perm(liste.rand, k_first=1)
per_k2 <- W_k_perm(liste.rand, k_first=2)
per_k3 <- W_k_perm(liste.rand, k_first=3)
per_k4 <- W_k_perm(liste.rand, k_first=4)
per_k5 <- W_k_perm(liste.rand, k_first=5)
per_k6 <- W_k_perm(liste.rand, k_first=6)
per_k7 <- W_k_perm(liste.rand, k_first=7)
per_k8 <- W_k_perm(liste.rand, k_first=8)
per_k9 <- W_k_perm(liste.rand, k_first=9)



# für die beobachten Daten:

W_k_obs <- function(obs.data, k_first){
  #~ auslagerbar
  is.wholenumber <- function(k_first, tol = .Machine$double.eps^0.5){ ### Abfangfunktion, ob k_first
    abs(k_first - round(k_first)) < tol}                                # ganzzahlig ist
  #~
  if(is.list(obs.data)!=TRUE){                                        ### Fehlerabfangung; Input muss eine 
                                                                        # Liste sein; jedes Element der Liste 
                                                                        # ist ein Dataframe/matrix mit SNPs
    stop("Stop! Input has to be a list. Each element hast to be a SNP-data-frame for a specific gen")
  }
  if(is.wholenumber(k_first)!=TRUE | k_first<=0){              ### Abfange; positiv und ganzzahlig?
    stop("Stop! k_first has to be a whole positiv number")
  }
  else{
    lapply(c(seq_along(obs.data)), 
           FUN=function(x){
             sorted.p <- sort(obs.data[[x]])[1:k_first]         ### Sortierung der p.-Werte und anschließend
             sorted.p <- sorted.p[!is.na(sorted.p)]   # nur die behalten, die auch wirkluch da sind (NAs raus)
             p.product <- prod(sorted.p)
           } #((())) Ende function in lapply
          ) #((())) Ende lapply
  } #((())) Ende esle
}  #((())) Ende Funktion
W_k_obs(liste.obs, 9)


### für beobachtete Daten; nur zum späteren Überprüfen
obs_k1 <- W_k_obs(liste.obs, k_first=1)
obs_k2 <- W_k_obs(liste.obs, k_first=2)
obs_k3 <- W_k_obs(liste.obs, k_first=3)
obs_k4 <- W_k_obs(liste.obs, k_first=4)
obs_k5 <- W_k_obs(liste.obs, k_first=5)
obs_k6 <- W_k_obs(liste.obs, k_first=6)
obs_k7 <- W_k_obs(liste.obs, k_first=7)
obs_k8 <- W_k_obs(liste.obs, k_first=8)
obs_k9 <- W_k_obs(liste.obs, k_first=9)

### k_per ist eine Liste mit den p-product-werten der der k-kleinsten snps 
  # für alle gene aus den permutierten daten
### k_obs ist eine Liste mit den p-product-werten der der k-kleinsten snps 
  # für alle gene
#################### WIRD DIESE FUNKTION EVTL NICHT MEHR GEBRAUCHT?
#################### Wenn doch, dann nur zur Überprüfung
p.W_k <- function(k_obs, k_per){
  lapply(c(seq_along(k_per)), 
         FUN=function(x){
           sum((k_per[[x]] <= k_obs[[x]]))/length(k_per[[x]])
         })
}


p.W_k(obs_k1, per_k1)
p.W_k(obs_k2, per_k2)
p.W_k(obs_k3, per_k3)
p.W_k(obs_k4, per_k4)
p.W_k(obs_k5, per_k5)


### mittels cbind wird eine matrix erstellt, die ZEILENWEISE die GENE enthält und 
  # als SPALTEN die UNTERSCHIEDLICHEN K
  # das hier ist zur Überprüfung, ob die Funktion funtioniert
cbind(unlist(p.W_k(obs_k1, per_k1)), unlist(p.W_k(obs_k2, per_k2)),  
      unlist(p.W_k(obs_k3, per_k3)), unlist(p.W_k(obs_k4, per_k4)),
      unlist(p.W_k(obs_k5, per_k5)), unlist(p.W_k(obs_k6, per_k6)),  
      unlist(p.W_k(obs_k7, per_k7)), unlist(p.W_k(obs_k8, per_k8)),
      unlist(p.W_k(obs_k9, per_k9)))





### Funktion um für verschiedene ks die p.Werte zu berechnen
p.W_ks <- function(ks, obs.data, perm.data){
        ### ks:         VEKTOR, der unterschiedliche Truncation Points (ks) enthält
        ### obs.data:   LISTE (!), die für jedes Gen/Set eine Element enthält, das
          #             wiederum als VEKTOR(!) die beobachteten (!) SNPs enthält
        ### perm.data:  LISTE (!), die für jedes Gen/Set eine Element enthält, das
          #             wiederum als MATRIX/DATAFRAME(!) die permutierten(!) SNPs enthält
  
  if(length(obs.data) != length(perm.data)){       ### Sind obs.data und perm.data gleich lang?
    stop("Number of variables/genes in observed and permutated data differ!")
  }
  
  else{
    
    if(is.null(nrow(perm.data[[1]]))){ ### wie viele Permutationen gab es; wird für p.W_ks gebraucht
      n.perm <- length(perm.data[[1]])   # es muss überprüft werden, ob das erste set nur aus einem
    }                                    # Vektor besteht. Wenn ja, muss dessen Laenge genommen werden
    else{
      n.perm <- nrow(perm.data[[1]])
    }
         
    products.perm <- lapply(ks, W_k_perm, perm.data=perm.data)     ### Funtion W_k_perm wird gebraucht !
                                                                     # Produkte der k-kleinsten SNPs
    products.obs <- lapply(ks, W_k_obs, obs.data=obs.data)         ### Funtion W_k_obs wird gebraucht !
                                                                     # Produkte der k-kleinsten SNPs
    
    n.gene <- length(obs.data)       ### Anzahl Gene muss bestimmt werden, damit für alle
                                        # ks über jedes Gen die Berechnung läuft
                                        
    
    temp <- lapply(seq_along(ks), FUN=function(x){
      unlist(lapply(c(1:n.gene), FUN=function(gene){   ### unlist, damit nicht listen in listen entstehen
        sum(products.perm[[x]][[gene]] <= products.obs[[x]][[gene]])/n.perm
      } #((())) Ende function vom inneren lapply
      ) #((())) Ende inneres lapply
      ) #((())) Ende unlist()
    } #((())) Ende function aeusseres lapply
    ) #((())) Ende auesseres lapply
    
    out <- matrix(unlist(temp), ncol=length(ks), byrow=FALSE)
    colnames(out) <- ks
    rownames(out) <- names(obs.data)
    
  } #((())) Ende else
  return(out)
} #((())) Ende Funktion


ps_per_k <- p.W_ks(c(1,5,8), obs.data=liste.obs, perm.data=liste.rand)


### Absichtlich zwei verschiedene Funktionen für p-werte über alle ks und nur die optimalen ks,
  # da so später ggf. noch mal geschaut werden kann, wie sich generell die 

### Funktion um die kleinsten ps für die ks zu finden

k.smallest_per_gen <- function(p.w_k_data){
        ### p.w_k_data: MATRIX, die die p-Werte der einzelnen ks enthält; Zeilen sind die Sets/Gene und
          #             Spalten sind die unterschiedlichen ks!
  
  if(is.matrix(p.w_k_data)!=TRUE){    ### p-Werte der unterschiedlichen ks müssen in Matrixform sein
    stop("Need matrix to proceed; rows have to be sets/genes, columns have to be different ks")
  }
  if(dim(p.w_k_data)[1] <= dim(p.w_k_data)[2]){   ### wenn es mehr Zeilen als Spalten gibt, koennte das der
                                                    # Hinweis sein, dass nicht die Zeilen die Gene und die
                                                    # Spalten die unterschiedlichen k sind; darum Warnung
    warning("There are more columns than rows. Are your rows the genes/sets and the columns the ks???")
  }
  
  smallest.p <- apply(p.w_k_data, MARGIN=1, min) ### für jedes Gen/Set den kleinsten p-Wert
  gene <- names(smallest.p)   ### Wie heißen die Gene/Sets?
  k <- sapply(seq_along(gene),  ### kurzes sapply, um ein VEKTOR zu erstellen, der den das
                                  # kleinste k für die kleinsten p-Werte eines Gens/Sets findet
              FUN=function(x){
                min(names(which((p.w_k_data[gene[x], ] == smallest.p[x])==TRUE))) ### name() damit der
              })                                          # Spaltenname(!) und nicht die Spaltenposition
                                                          # als Bestimmung von k ausgewählt wird
  out <- data.frame(gene, k)
  return(out)
}

selected_ks <- k.smallest_per_gen(ps_per_k)



### auswahl der entsprechenden snps in der Liste mit den p-Werten der
  # beobachteten SNPs

selected.snps.obs <- function(obs.data, selected_ks){
  out <- lapply(seq_along(obs.data), 
         FUN=function(var){
           sort(obs.data[[selected_ks[var,1]]])[1:selected_ks[var,2]]  
         })
  names(out) <- names(obs.data)
  return(out)
}

OBS.snps <- selected.snps.obs(obs.data=liste.obs, selected_ks=selected_ks)



selected.snps.perm <- function(perm.data, selected.snps.obs){
  out <- lapply(seq_along(perm.data),
                FUN=function(var){
                   if(length(dim(perm.data[[var]]))< 2){
                     temp.data <- as.matrix(perm.data[[var]], ncol=1)
                     colnames(temp.data) <- names(selected.snps.obs[[var]])
                     temp.data
                   }else{
                     temp.data <- perm.data[[var]][,which(colnames(perm.data[[var]]) %in%  names(selected.snps.obs[[var]]))]
                     if(length(dim(temp.data)) <2){
                       temp.data <- as.matrix(temp.data, ncol=1)
                       colnames(temp.data) <- names(selected.snps.obs[[var]])
                     }
                    
                     temp.data
                   }
                     
                })
  names(out) <- names(selected.snps.obs)
  return(out)
}

PERM.snps <- selected.snps.perm(perm.data=liste.rand, selected.snps.obs=OBS.snps)


colnames(PERM.snps[[1]])==sort(names(OBS.snps[[1]]))
colnames(PERM.snps[[2]])==sort(names(OBS.snps[[2]]))
colnames(PERM.snps[[3]])==sort(names(OBS.snps[[3]]))
colnames(PERM.snps[[4]])==sort(names(OBS.snps[[4]]))

table(PERM.snps[[1]]==liste.rand[[1]][,colnames(PERM.snps[[1]])])
table(PERM.snps[[2]]==liste.rand[[2]][,colnames(PERM.snps[[2]])])
table(PERM.snps[[3]]==liste.rand[[3]][,colnames(PERM.snps[[3]])])


####
# to do:
### simulation für gsea-snp
### wiesowählen bei der ssea-methode die autoren erst alle snps aus und nicht nur die, die auch für die eigentlichen PW interessant
  # sind? shcließlich ist das ja nur für die wahl der optimalen snps und nicht als reine ARTP gedacht; Ronja fragen!



