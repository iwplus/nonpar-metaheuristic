##############################################################
###### Library/Packages yang diperlukan untuk run program ###########
##############################################################

library(pracma)

##############################################################
############ Fungsi tambahan #################################
##############################################################



############## Fungsi menghapus unsur nol di vektor #######

removezero = function(a)
{
 indzero = 0
 ind = 1
 while (ind <= length(a)) #### cari indeks entri tak-nol pertama
 {
	if (a[ind] != 0) 
	{
		indzero = ind
		ind = length(a)+5
	}
 ind = ind + 1
 }
 
 for (i in (indzero+1):length(a)) ### cari semua indeks entri tak-nol
 {
	if (a[i] != 0)
	{
		indzero = cbind(indzero,i)
	}
 }
 
 b = matrix(0,nrow=1,ncol=length(indzero))
 
 for (j in 1:length(b))
 {
	b[j] = a[indzero[j]]
 }
 
 return(b)
}

###################################################################################################


###########################################
### (1) Bagian Genetics Algorithm ####
###########################################

###############################################
####Fungsi-fungsi yang berguna #########
###############################################



### fungsi seleksi populasi menggunakan rank selection

pop_selection <- function(fitness,pop,n_select){
    #print(fitness)
    sortfitness = sort(fitness, index.return = TRUE);
    len_fitness = nrow(pop)
    prob_fitness = matrix(0, nrow=1, ncol=len_fitness)

    # print("panjang fitness")
    # print(len_fitness)
    # print(nrow(pop))
    # print("fitness")
    # print(fitness)

    rank_sum <- len_fitness*(len_fitness+1)/2
    for (i in 1:len_fitness){
        prob_fitness[1,sortfitness$ix[i]] <- i/rank_sum
    }
    # print("prob fitness")
    # print(prob_fitness)

    # prob_fitness = matrix(0, nrow=1, ncol=len_fitness)
    # for (i in 1:len_fitness){
    #     prob_fitness[1,i] <- rank_weight[1,i]/rank_sum
    # }

    cum_prob = matrix(0,nrow=1,ncol=len_fitness)

    for (i in 1:len_fitness){
        if (i==1){
            cum_prob[1,i] <- prob_fitness[1,i]
        }
        else{
            cum_prob[1,i] <- cum_prob[1,i-1] + prob_fitness[1,i]
        }
    }
    # print("cum_prob : ")
    # print(cum_prob)

    pop_selected = matrix(0, ncol=ncol(pop),nrow=n_select)

    for (i in 1:n_select){
        p = runif(1)
        for (j in 1:len_fitness){
            if (cum_prob[1,j]>= p){
                ind_selected <- j
                #print(ind_selected)
                break
            }
        }
        pop_selected[i,] <- pop[ind_selected,]
    }
    return(pop_selected)
}

## fungsi cross-over populasi

swap_suffix <- function(a,b){
    s = sample(1:length(a),1)
    off_spr1 <- a
    off_spr2 <- b

    off_spr1[s:length(a)] = b[s:length(a)]
    off_spr2[s:length(a)] = a[s:length(a)]

    off_spr = rbind(off_spr1,off_spr2)
    return(off_spr)
}

cross_over <- function(p_cross,nvarkernel,nvarspline,nknot,pop_selected){
    n_selected = nrow(pop_selected)
    n_col = nvarkernel + (nknot*nvarspline)
    pop_cross = matrix(0, ncol=n_col,nrow=2*n_selected)
    for (i in 1:n_selected){
        p = runif(1)
        if (p_cross > p){
            if (i == n_selected){
                ind_1 = i
                ind_2 = 1
            }
            else{
                ind_1 = i
                ind_2 = i+1
            }
            a_1 = pop_selected[ind_1,]
            a_2 = pop_selected[ind_2,]

            ## cross-over pada bandwidth
            band_cross = swap_suffix(a_1[1:nvarkernel],a_2[1:nvarkernel])

            ## cross-over pada knot
            knot_cross = swap_suffix(a_1[(nvarkernel+1):n_col],a_2[(nvarkernel+1):n_col])

            ## satukan kedua hasil cross-over
            off_spr = cbind(band_cross,knot_cross)

            ## simpan ke dalam populasi hasil cross-over
            pop_cross[((2*(i-1))+1),] = off_spr[1,]
            pop_cross[(2*i),] = off_spr[2,]


        }
    }

    return(pop_cross)

}

## Fungsi mutasi pada off-spring

mutasi <- function(p_mutasi,pop_cross,batasvarkernel,batasvarspline,nknot){
    for (i in 1:nrow(pop_cross)){
        p = runif(1)
        if (p_mutasi > p){
            ## ambil secara random indeks untuk mutasi bandwidth dan knot
            #print(nrow(batasvarkernel))
            if (is.null(nrow(batasvarkernel))){
                ind_b <- 1
            }
            else{
                ind_b <- sample(1:nrow(batasvarkernel))
            }
            
            if (is.null(nrow(batasvarspline))){
                ind_s = 1
            }
            else{
                ind_s <- sample(1:nrow(batasvarspline))
            }
            
            ## mutasi bandwidth
            #if (nrow(batasvarkernel) == 1){
            if (is.null(nrow(batasvarkernel))){
                new_band = runif(1, min=batasvarkernel[1], max=batasvarkernel[2])
            }
            else{
                new_band = runif(1, min=batasvarkernel[ind_b,1], max=batasvarkernel[ind_b,2])
            }
            pop_cross[i,ind_b] = new_band

            ## mutasi knot
            for (j in 1:nknot){
                #if (nrow(batasvarspline) == 1){
                if (is.null(nrow(batasvarspline))){
                    new_knot = runif(1, min=batasvarspline[1], max=batasvarspline[2])
                }
                else{
                    new_knot = runif(1, min=batasvarspline[ind_s,1], max=batasvarspline[ind_s,2])
                }
                pop_cross[i,nrow(batasvarkernel)+ind_s+(j-1)] = new_knot
            }
        }
    }
return(pop_cross)
}

## Fungsi inisialisasi populasi

inisialisasi_pop <- function(N,nvarkernel,batasvarkernel,nvarspline,nknot,batasvarspline){
    X1 = matrix(0, nrow=N, ncol = nvarkernel); ################ bangkitkan bandwidth awal
    for (i in 1:N)
    {
    for (j in 1:nvarkernel)
    {
        #indrn = sample.int(100,1)
        #print(batasvarkernel[j,1])
        if (nvarkernel == 1){
            #X1[i,j] = batasvarkernel[1]+((indrn-1)*(batasvarkernel[2]-batasvarkernel[1])/npartisidomainvar)
            X1[i,j] = runif(1, min=batasvarkernel[1], max=batasvarkernel[2])
        } 
        else{
            #X1[i,j] = batasvarkernel[j,1]+((indrn-1)*(batasvarkernel[j,2]-batasvarkernel[j,1])/npartisidomainvar)
            X1[i,j] = runif(1, min=batasvarkernel[j,1], max=batasvarkernel[j,2])
        }
            
    }
    }

    ###
    nvarspline = ncol(t)
    batasvarspline = batasvar[(nvarkernel+1):nrow(batasvar),]
    X2 = matrix(0, nrow=N, ncol = nvarspline*nknot); ################# bangkitkan knot awal
    for (i in 1:N)
    {
    for (j in 1:nvarspline)
    {
        #rn = sample.int(100,nknot)
        for (k in 1:nknot)
        {
            if (nvarspline == 1){
                #X2[i,((j-1)*nknot)+k] = batasvarspline[1]+((rn[k]-1)*(batasvarspline[2]-batasvarspline[1])/npartisidomainvar)
                X2[i,((j-1)*nknot)+k] = runif(1, min=batasvarspline[1], max=batasvarspline[2])
            }
            else{
                #X2[i,((j-1)*nknot)+k] = batasvarspline[j,1]+((rn[k]-1)*(batasvarspline[j,2]-batasvarspline[j,1])/npartisidomainvar)
                X2[i,((j-1)*nknot)+k] = runif(1, min=batasvarspline[j,1], max=batasvarspline[j,2])
            }	
        }
            
            
    }
    }

    X = cbind(X1,X2) #################### Augment bandwidth dan knot

    return(X)
}

#########################################
#### Inisialisasi input GA* ############
#########################################
### *sebagai contoh saja (di masing-masing bagian sudah ditaruh inisialisasi juga) ################

#N = 5; ### Banyaknya calon solusi awal
#itermaks = 100; ### banyaknya iterasi maksimum
#nvar = 5; ### diisi dengan banyaknya variabel 


################################################
############ Fungsi AG ########################
################################################

 
genetics_alg <- function(N,itermax,f,y,x,t,batasvar,degree,nknot) #### bagian awal fungsi pso yang meminimumkan f
{
 ## parameter GA lainnya
 p_cross = 0.7
 p_mutasi = 0.6
 
 #npartisidomainvar = 100 #### banyaknya partisi sama besar pada domain variabel
 n_select = floor(N/3)
 
 ## banyaknya var independent dan batasnya untuk regresi kernel
 nvarkernel = ncol(x)
 batasvarkernel = batasvar[1:nvarkernel,]
 #print(batasvarkernel[1])
 
 ## banyaknya var independent dan batasnya untuk regresi spline truncated
 nvarspline = ncol(t)
 batasvarspline = batasvar[(nvarkernel+1):nrow(batasvar),]

 ## Inisialisasi populasi
 pop_baru = inisialisasi_pop(N,nvarkernel,batasvarkernel,nvarspline,nknot,batasvarspline)
 #print(X)
  
 for (iter in 1:itermaks){
    pop <- pop_baru
    
    # print("banyak populasi")
    # print(nrow(pop))


    fitness = matrix(0,nrow=1,ncol=nrow(pop));
    for (i in 1:nrow(pop)){
        Xtemp <- pop[i,]
        val <- f(y,x,t,Xtemp,degree)
        if (is.nan(val)){
            fitness[1,i] <- 0;
        }
        else{
            fitness[1,i] <- 1/(val+1);
        }
        
    } 
    #print(fitness)    
    ## seleksi populasi
    pop_selected <- pop_selection(fitness,pop,n_select)

    ## cross-over populasi hasil seleksi
    pop_cross <- cross_over(p_cross,nvarkernel,nvarspline,nknot,pop_selected)  ## note: jangan lupa sediakan nilai p_cross

    ## mutasi off-spring hasil cross-over
    pop_mutasi <- mutasi(p_mutasi,pop_cross,batasvarkernel,batasvarspline,nknot) ## note: jangan lupa sediakan nilai p_mutasi

    ## gabungkan off-spring dengan populasi yang terseleksi
    pop_baru <- rbind(pop_selected,pop_mutasi)  ## perlu dicek apakah ini menyebabkan error

 }
 fitness = matrix(0,nrow=1,ncol=nrow(pop_baru));
 for (i in 1:nrow(pop_baru)){
    Xtemp = pop_baru[i,]
    val <- f(y,x,t,Xtemp,degree)
    if (is.nan(val)){
        fitness[1,i] <- 0;
    }
    else{
        fitness[1,i] <- 1/(val+1);
    }
 }
 sortfitness = sort(fitness, index.return = TRUE);
 best_indeks = sortfitness$ix[nrow(pop_baru)]
 g <- pop_baru[best_indeks,]
#  c(g,f(y,x,t,g,degree))
 c(g,fitness[1,best_indeks])

} ######################################### akhir dari fungsi genetics algorithm

#output = genetics_alg(N,itermaks,gcvmix,y,x,t,batasvar,degree,nknot) ##################### Contoh cara panggil fungsi optimisasi GA

#########################################################
#########################################################


#####################################################################################
############# (2) Bagian Regresi Kernel + Spline Truncated ###############################################
#####################################################################################

################################# Fungsi-fungsi penting #####################
#############################################################################

#### Fungsi kernel ############### 

kernel = function(u) #### kernel bi-square
{
 if (abs(u)<=1) {ker = 0.9375*(1-u^2)^2}
 else {ker = 0}
 return(ker)
}

########## Norma Euclid ##############

euclidnorm = function(b) #### Norma Euclid 
{
 norma = 0
 for (i in 1:length(b))
 {
	norma = norma + b[i]^2
 }
 return(norma)
}#######################################


########### Fungsi W pada estimator Nadaraya-Watson ################

nadwatson = function(x0,x1,x,h) ###### Fungsi W pada estimator Nadaraya-Watson
{
 m = length(x)
 K = 0

 for (s in 1:m)
 {
	if (h == 0){
        Ktemp = 0
        }
	else {Ktemp = kernel(((x0-x[s])/h))}
	K = K + Ktemp
 }
 
 if (h == 0){
    Ktemp2 = 0
    }
 else {
    Ktemp2 = kernel(((x0-x1)/h))
    }
 
 W = Ktemp2/K
 return(W)
}################################################################################

###################### Fungsi untuk menghitung matriks D(alpha) ##############

dalpha = function(x,h)
{
	n = nrow(x)
	Dalpha = matrix(0, nrow = n, ncol = n)

	#for (indvar in 1:nrow(indeks))
	for (indvar in 1:nrow(indeksx))
	{
 		D = matrix(0, ncol = n, nrow = n) ### Inisialisasi matriks D(h)
 
 		for (k in 1:n) ##### pengisian entri matriks D(h) dengan fungsi W dari estimator Nadaraya-Watson
 		{
   			for (l in 1:n)
   			{
				D[k,l] = 1/n*(nadwatson(x[k,indvar],x[l,indvar],x[,indvar],h[indvar]))
   			}
 		}
 		Dalpha = Dalpha + D
	}
	return(Dalpha)
}

############ Fungsi truncated ############

truncate = function(t,k,d)
{
  if (t >= k)
   {
   	trun = (t-k)^d 
   }
  else
   {
	trun = 0
   }
 return(trun)
}

################################################################################


################################################################################
################### Fungsi GCV untuk Mixed ######################################
################################################################################


gcvmix <- function(y,x,t,hknot,degree) #### hknot = [bandwidth | knot] (augmented)
{
 ################## Bagian Regresi Kernel #####################
 
 nvarkernel = ncol(x)
 h = hknot[1:nvarkernel] ##### Bandwidth
 Dalpha = dalpha(x,h) ###### hitung matriks D(alpha)
 
 ##############################################################

 ############### Bagian Spline Truncated ######################

 knot = hknot[(nvarkernel+1):length(hknot)] ##### Knot
 #print(knot)
 #print(length(knot))
 n = nrow(t)
 nmat = ncol(t)
 #print(nmat)
 jumknot = length(knot)/nmat
 #print(nknot)
 m = degree + 1 + jumknot
 Gk = matrix(0,nrow = n, ncol = m) ###### inisialisasi matriks G(k)
 
 for (indcol in 1:m) ###### menghitung G(k_1)
 {
	if (indcol <= degree+1)
	{
	  	for (indrow in 1:n)
	  	{
			Gk[indrow,indcol] = t[indrow,1]^(indcol-1)
	  	}
	}
	else
	{
		for (indrow in 1:n)
		{
			Gk[indrow,indcol] = truncate(t[indrow,1],knot[indcol-(degree+1)],degree)
		}	
	}
 }	
 
 #print(Gk) ###### coba print
if (nmat >= 2){
 for (indvar in 2:nmat) #### menghitung G(k_2),...,G(k_p) dan digabung jadi G(k)
 {
	Gtemp = matrix(0, nrow = n, ncol = m)
	for (indcol in 1:m) ###### menghitung G(k_1)
 	{
		if (indcol <= (degree+1))
		{
	  		for (indrow in 1:n)
	  		{
				
				Gtemp[indrow,indcol] = t[indrow,indvar]^(indcol-1)
	  		}
		}
		else
		{
			for (indrow in 1:n)
			{
				Gtemp[indrow,indcol] = truncate(t[indrow,indvar],knot[((indvar-1)*jumknot)+(indcol-(degree+1))],degree)
			}	
		}
 	}	

 	Gk = cbind(Gk,Gtemp)	
 }
}
 #print(Gk)
 #print(t(Gk)%*%Gk)
 I = diag(nrow(Dalpha))
 K = Gk%*%pinv((t(Gk)%*%Gk))%*%t(Gk)%*%(I-Dalpha)
 M = K + Dalpha
 H = I-M
 HA = 1/n*(euclidnorm(H%*%y))^2
 HB = (1/n*sum(diag(H)))^2
 gcvm = HA/HB
 #print(gcvm)
 return(gcvm)
}


################################################################
####### (3) Pengolahan data (eksekusi) #############
################################################################

##### Baca data dan parameter lain yang diperlukan #############

mydata = read.csv('/Users/irwansyah/Documents/mixed-spline-project/experiment/a_data.csv',header=TRUE)
dataindeks = read.csv('/Users/irwansyah/Documents/mixed-spline-project/experiment/indeks_a.csv',header=TRUE)


####### Data (variabel) respon/dependent ###########

var_ind_start = 3 ## input indeks kolom awal variabel independen - 1
var_dep_start = 3 ## input indeks kolom awal variabel dependen

mydata = data.frame(mydata)
y = as.matrix(mydata[,var_dep_start]) 

### Input indeks variabel yang sudah tercatat di file csv indeks####

dataindeks = data.frame(dataindeks)
indeks = as.matrix(dataindeks)
print(indeks[,1])

##### Data variabel x ########################

indeksx = removezero(indeks[,1]) ##### indeks variabel x ada di kolom ke-1 di file cvs indeks 
print(length(indeksx))
x = matrix(ncol=length(indeksx),nrow=nrow(y))

for (i in 1:length(indeksx))
{
 for (j in 1:nrow(y))
 {
   x[j,i] = mydata[j,indeksx[i]+var_ind_start] #### variabel respon dimulai dari kolom ke-var_ind_start+1 di file csv data
 }
}
print(x)
print(ncol(x))


###### Ambil data variabel t ################

indekst = removezero(indeks[,2]) ##### indeks variabel x ada di kolom ke-2 di file csv indeks

t = matrix(ncol=length(indekst),nrow=nrow(y))

for (i in 1:length(indekst))
{
 for (j in 1:nrow(y))
 {
   t[j,i] = mydata[j,indekst[i]+var_ind_start] #### variabel respon dimulai dari kolom ke-var_ind_start + 1 di file csv data
 }
}
# print(t)
# print(ncol(t))
###### batas atas dan bawah masing-masing bandwidth ########

 batasvarkernel = matrix(0, nrow = ncol(x), ncol = 2)

 for (i in 1:ncol(x))
 {
	batasvarkernel[i,1] = 0.1
	batasvarkernel[i,2] = 100
 }

#print(batasvarkernel)
############ batas atas dan bawah masing-masing knot ########

batasvarspline = matrix(0,nrow = ncol(t), ncol = 2)

for (i in 1:ncol(t))
{
	batasvarspline[i,1] = min(t[,i])
	batasvarspline[i,2] = max(t[,i])
}
#print(batasvarspline)
#### Inisialisasi input PSO untuk Regresi Kernel##############

batasvar = rbind(batasvarkernel,batasvarspline)
#print(batasvar)

N = 20; ### Banyaknya calon solusi awal
itermaks = 10; ### banyaknya iterasi maksimum

degree = 2 ###### derajat polinom spline
nknot = 3 ######## banyaknya knot

################################################################

bandknotmgcv = genetics_alg(N,itermaks,gcvmix,y,x,t,batasvar,degree,nknot); ######## pencarian bandwidth optimum dengan GA 

bandknotmgcv

########################## Cetak matriks D(alpha), M(k), dan data hasil prediksi #################

bandakhir = bandknotmgcv[1:ncol(x)]
Dalphaf = dalpha(x,bandakhir)

knotf = bandknotmgcv[(ncol(x)+1):(length(bandknotmgcv)-1)] ##### Knot
 n = nrow(t)
 nmat = ncol(t)
 jumknot = length(knotf)/nmat
 m = degree + 1 + jumknot
 Gk = matrix(0,nrow = n, ncol = m) ###### inisialisasi matriks G(k)
 
 for (indcol in 1:m) ###### menghitung G(k_1)
 {
	if (indcol <= degree+1)
	{
	  	for (indrow in 1:n)
	  	{
			Gk[indrow,indcol] = t[indrow,1]^(indcol-1)
	  	}
	}
	else
	{
		for (indrow in 1:n)
		{
			Gk[indrow,indcol] = truncate(t[indrow,1],knotf[indcol-(degree+1)],degree)
		}	
	}
 }	
 
 #print(Gk) ###### coba print
if (nmat>=2){
 for (indvar in 2:nmat) #### menghitung G(k_2),...,G(k_p) dan digabung jadi G(k)
 {
	Gtemp = matrix(0, nrow = n, ncol = m)
	for (indcol in 1:m) ###### menghitung G(k_1)
 	{
		if (indcol <= (degree+1))
		{
	  		for (indrow in 1:n)
	  		{
				
				Gtemp[indrow,indcol] = t[indrow,indvar]^(indcol-1)
	  		}
		}
		else
		{
			for (indrow in 1:n)
			{
				Gtemp[indrow,indcol] = truncate(t[indrow,indvar],knotf[((indvar-1)*jumknot)+(indcol-(degree+1))],degree)
			}	
		}
 	}	

 	Gk = cbind(Gk,Gtemp)	
 }
}
 #print(Gk)
 #print(t(Gk)%*%Gk)
 I = diag(nrow(Dalphaf))
 K = Gk%*%pinv((t(Gk)%*%Gk))%*%t(Gk)%*%(I-Dalphaf)
 M = K + Dalphaf

ypred = M%*%y
print(ypred)

#plot(y,type="b", col='red'); plot(ypred, col = 'blue', add = TRUE)
#plot(y,type="b",col="red"); lines(ypred,col="green")
output = cbind(y,ypred); 
write.csv(output, "/Users/irwansyah/Documents/mixed-spline-project/experiment/result_a_ga_x1-5-kernel_x6-9-spline.csv")

parameters = c(degree,nknot,knotf,bandakhir);
write.csv(parameters, "/Users/irwansyah/Documents/mixed-spline-project/experiment/params_a_ga_x1-5-kernel_x6-9-spline.csv")
# library(ggplot2)
# ggplot(ypred)

### hitung nilai GCV mix optimal

I = diag(nrow(Dalphaf))
H = I-M
HA = 1/n*(euclidnorm(H%*%y))^2
HB = (1/n*sum(diag(H)))^2
gcvm = HA/HB
print("GCV mix : ")
print(gcvm)

## track the GCV mix values
