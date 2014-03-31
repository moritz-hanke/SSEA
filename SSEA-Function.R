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





k_counter <- NULL

sortiererei <- function(daten, k_first){
  #~ auslagerbar
  is.wholenumber <- function(k_first, tol = .Machine$double.eps^0.5){ ### Abfangfunktion, ob k_first
    abs(k_first - round(k_first)) < tol}                                # ganzzahlig ist
  #~
  if(is.list(daten)!=TRUE){                                           ### Fehlerabfangung; Input muss eine Liste sein
    print("Stop! Input has to be a list. Each element hast to be a SNP-data-frame for a specific gen")
  }
  if(is.wholenumber(k_first)!=TRUE | k_first<=0){
    print("Stop! k_first has to be a whole positiv number")
  }
  else{
    lapply(c(seq_along(daten)), FUN=function(var){   ### äußere lapply für Genauswahl der in der Liste
      apply(daten[[var]], MARGIN=1,             ### inneres apply für zeilenweisen Zugriff auf die 
                                                       # Datensätze/matrizen mit den Variablen
            #~ auslagerbar
            FUN=function(x){
              sorted.p <- sort(x)[1:k_first]         ### Sortierung der p.-Werte und anschließend
              sorted.p <- sorted.p[!is.na(sorted.p)]   # nur die behalten, die auch wirkluch da sind (NAs raus)
              #k_counter
              p.product <- prod(sorted.p)
              
              
            } ### Ende der Funktion vom inneren apply
            #~
      ) ### Ende inneres apply
    } ### Ende function von äußerem lapply
    ) ### Ende äußeres lapply
  } ### Ende von else
  
} ## Ende der Funktion

d <- sortiererei(liste.rand, 9)

d[[3]][100]
d

# für permutierte dateb
per_k1 <- sortiererei(liste.rand, k_first=1)
per_k2 <- sortiererei(liste.rand, k_first=2)
per_k3 <- sortiererei(liste.rand, k_first=3)
per_k4 <- sortiererei(liste.rand, k_first=4)
per_k5 <- sortiererei(liste.rand, k_first=5)




# für die beobachten Daten, MUSS noch berarbeitet werden !!!:

sortierung.obs <- function(daten.obs, k_first){
  lapply(c(seq_along(daten.obs)), 
         FUN=function(x){
           sorted.p <- sort(daten.obs[[x]])[1:k_first]         ### Sortierung der p.-Werte und anschließend
           sorted.p <- sorted.p[!is.na(sorted.p)]   # nur die behalten, die auch wirkluch da sind (NAs raus)
           #k_counter
           p.product <- prod(sorted.p)
})}
sortierung.obs(liste.obs, 9)

# für permutierte dateb
obs_k1 <- sortierung.obs(liste.obs, k_first=1)
obs_k2 <- sortierung.obs(liste.obs, k_first=2)
obs_k3 <- sortierung.obs(liste.obs, k_first=3)
obs_k4 <- sortierung.obs(liste.obs, k_first=4)
obs_k5 <- sortierung.obs(liste.obs, k_first=5)



sum((d_k1[[1]] <= e_k1[[1]]))/length(d_k1[[1]])
sum((d_k2[[1]] <= e_k2[[1]]))/length(d_k2[[1]])
sum((d_k3[[1]] <= e_k3[[1]]))/length(d_k3[[1]])
sum((d_k4[[1]] <= e_k4[[1]]))/length(d_k4[[1]])
sum((d_k5[[1]] <= e_k5[[1]]))/length(d_k5[[1]])




### k_per ist eine Liste mit den p-product-werten der der k-kleinsten snps 
  # für alle gene aus den permutierten daten
### k_obs ist eine Liste mit den p-product-werten der der k-kleinsten snps 
  # für alle gene
p.W_k <- function(k_obs, k_per){
  lapply(c(seq_along(k_per)), 
         FUN=function(x){
           sum((k_per[[x]] <= k_obs[[x]]))/length(k_per[[x]])
         })
}


unlist(p.W_k(obs_k1, per_k1))



lapply(c(seq_along(d_k1)), 
       FUN=function(x){
         sum((d_k1[[x]] <= e_k1[[x]]))/length(d_k1[[x]])
       })




