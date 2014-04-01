# Funktionen, um SSEA durchzuführen


p.rand <- matrix(round(runif(2000, 0, 1), 8), ncol=20, byrow=T)
p.obs <- matrix(round(runif(20, 0, 1), 8), ncol=20, byrow=T)
colnames(p.rand) <- letters[1:20]
colnames(p.obs) <- letters[1:20]


set1 <- letters[c(1:4, 6)]
set2 <- letters[c(5,7:13)]
set3 <- letters[14:20]


liste.obs <- vector(mode="list", length=3)
liste.obs[[1]] <- p.obs[,colnames(p.obs) %in% set1]
liste.obs[[2]] <- p.obs[,colnames(p.obs) %in% set2]
liste.obs[[3]] <- p.obs[,colnames(p.obs) %in% set3]
names(liste.obs) <- c("set1", "set2", "set3")


liste.rand <- vector(mode="list", length=3)
liste.rand[[1]] <- p.rand[,colnames(p.rand) %in% set1]
liste.rand[[2]] <- p.rand[,colnames(p.rand) %in% set2]
liste.rand[[3]] <- p.rand[,colnames(p.rand) %in% set3]
names(liste.rand) <- c("set1", "set2", "set3")

#################
# BEISPIEL für Schachtelung von applys u.Ä.
# für die zufälligen/permutierten Daten soll eine Sortierung für
# jedes set/gen der variablen/snps durchgenommen werden
b <- lapply(c(1,2,3), FUN=function(var){    # äußere lapply für Genauswahl der in der Liste
  apply(liste.rand[[var]], MARGIN=1,        # inneres apply für zeilenweisen Zugriff auf die 
                                            # Datensätze/matrizen mit den Variablen
      FUN=function(x){
        sort(x)[]
      }
    )
  }
)





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
      apply(perm.data[[var]], MARGIN=1,             ### inneres apply für zeilenweisen Zugriff auf die 
                                                       # Datensätze/matrizen mit den Variablen
            #~ auslagerbar
            FUN=function(x){
              sorted.p <- sort(x)[1:k_first]         ### Sortierung der p.-Werte und anschließend
              sorted.p <- sorted.p[!is.na(sorted.p)]   # nur die behalten, die auch wirkluch da sind (NAs raus)
              p.product <- prod(sorted.p)
              
              
            } ### Ende der Funktion vom inneren apply
            #~
      ) ### Ende inneres apply
    } ### Ende function von äußerem lapply
    ) ### Ende äußeres lapply
  } ### Ende von else
  
} ## Ende der Funktion

d <- W_k_perm(liste.rand, 9)

d[[3]][100]
d

# für permutierte daten; nur zum späteren Überprüfen
per_k1 <- W_k_perm(liste.rand, k_first=1)
per_k2 <- W_k_perm(liste.rand, k_first=2)
per_k3 <- W_k_perm(liste.rand, k_first=3)
per_k4 <- W_k_perm(liste.rand, k_first=4)
per_k5 <- W_k_perm(liste.rand, k_first=5)




# für die beobachten Daten, MUSS noch berarbeitet werden !!!:

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
           })
  }
}
W_k_obs(liste.obs, 9)

### für beobachtete Daten; nur zum späteren Überprüfen
obs_k1 <- W_k_obs(liste.obs, k_first=1)
obs_k2 <- W_k_obs(liste.obs, k_first=2)
obs_k3 <- W_k_obs(liste.obs, k_first=3)
obs_k4 <- W_k_obs(liste.obs, k_first=4)
obs_k5 <- W_k_obs(liste.obs, k_first=5)


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
cbind(unlist(p.W_k(obs_k1, per_k1)), unlist(p.W_k(obs_k2, per_k2)))


#es wird eine funktion gebraucht, die die unterschiedlichen p.W_k auf einmal erzeugt
lapply(c(1,2), W_k_obs, daten.obs=liste.obs)       ### Zeilen sind Gene, Spalten unterschiedliche ks
cbind(unlist(obs_k1),unlist(obs_k2),unlist(obs_k3))


products.perm <- lapply(c(1:5), W_k_perm, daten=liste.rand)
products.obs <- lapply(c(1:5), W_k_obs, daten.obs=liste.obs)


n.perm <-100

# Funktion um für verschiedene ks die p.Werte zu berechnen
p.W_ks <- function(ks, obs.data, perm.data){
  
  if(length(obs.data) != length(perm.data)){
    stop("Number of variables/genes in observed and permutated data differ!")
  }
  
  # Abfangen von Fehlern muss hier her!!!!
  
  products.perm <- lapply(ks, W_k_perm, perm.data=perm.data)     ### Funtion W_k_perm wird gebraucht !
                                                                # Produkte der k-kleinsten SNPs
  products.obs <- lapply(ks, W_k_obs, obs.data=obs.data)    ### Funtion W_k_obs wird gebraucht !
                                                                # Produkte der k-kleinsten SNPs
  
  n.gene <- length(liste.obs)       ### Anzahl Gene muss bestimmt werden, damit für alle
                                      # ks über jedes Gen die berechnung läuft
  
  lapply(ks, FUN=function(x){
    unlist(lapply(c(1:n.gene), FUN=function(gene){
      sum(products.perm[[x]][[gene]] <= products.obs[[x]][[gene]])/n.perm
      }
      )
    )  
  }
  )
}

p.W_ks(c(1,2,3))




str(test[[1]])


products.perm[[2]][[3]] <= products.obs[[2]][[3]]


products.perm[[1]] <= products.obs


sum((lapply(c(1,2), W_k_perm, daten=liste.rand)[[2]][[3]] <= lapply(c(1,2), W_k_obs, daten.obs=liste.obs)[[2]][[3]]))/
  nrow(liste.rand[[1]])









lapply(c(seq_along(d_k1)), 
       FUN=function(x){
         sum((d_k1[[x]] <= e_k1[[x]]))/length(d_k1[[x]])
       })





