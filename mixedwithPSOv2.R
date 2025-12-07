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
### (1) Bagian particle swarm optimization ####
###########################################

###############################################
####Daftar fungsi-fungsi yang berguna #########
###############################################


### fungsi bobot inersia ####

omega <- function(t,tmaks){
	bobinersia = 0.9-(((0.9-0.4)/tmaks)*t);
	bobinersia
}

### fungsi pembangkit matriks R1 dan R2 ##########

generater12 <- function(nvar){
 entridiag = c(matrix(runif(1*nvar), ncol = nvar));
 if (nvar == 1)
 {
	matriksr = entridiag;
 }
 else
 {
	matriksr = diag(entridiag);
 }
 return(matriksr)
}



#########################################
#### Inisialisasi input PSO* ############
#########################################
### *sebagai contoh saja (di masing-masing bagian sudah ditaruh inisialisasi juga) ################

#N = 5; ### Banyaknya calon solusi awal
#itermaks = 100; ### banyaknya iterasi maksimum
#nvar = 5; ### diisi dengan banyaknya variabel 


################################################
############ Fungsi PSO ########################
################################################

 
pso <- function(N,itermax,f,y,x,t,batasvar,degree,nknot) #### bagian awal fungsi pso yang meminimumkan f
{
 npartisidomainvar = 100 #### banyaknya partisi sama besar pada domain variabel
 
 nvarkernel = ncol(x)
 batasvarkernel = batasvar[1:nvarkernel,]
 #print(batasvarkernel[1])
 X1 = matrix(0, nrow=N, ncol = nvarkernel); ################ bangkitkan bandwidth awal
 for (i in 1:N)
 {
   for (j in 1:nvarkernel)
   {
	indrn = sample.int(100,1)
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
	rn = sample.int(100,nknot)
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

 #print(X)
 nvarknot = nvarkernel + (nvarspline*nknot) 
 p = X; ### personal best awal
 v = matrix(0,N,nvarknot); ### kecepatan awal


 fitness = matrix(0,1,N);

 
 for (i in 1:N){
	Xtemp = X[i,]
	fitness[i] = 1/(f(y,x,t,Xtemp,degree)+1);
 } 

 sortfitness = sort(fitness, index.return = TRUE);

 g = X[sortfitness$ix[1],]; ### solusi minimum

 #print(fitness)
 #print(X)

 for (iter in 2:itermaks){


 ### update kecepatan particle ###

 for (i in 1:N){
 	v[i,] = (omega(iter,itermaks)*v[i,])+(2*(p[i,]-X[i,])%*%generater12(nvarknot))+(2*(g-X[i,])%*%generater12(nvarknot));
 }


 ##### update solusi #####

 for (i in 1:N){
	X[i,] = X[i,] + v[i,];
	Xtemp2 = X[i,]
	if (f(y,x,t,Xtemp2,degree)< fitness[i]){
		p[i,] = X[i,];
	}
	
	if (f(y,x,t,X[i,],degree)< f(y,x,t,g,degree)){
		g = X[i,];
	}
 }

 for (i in 1:N){
	Xtemp3 = X[i,]
	fitness[i] = f(y,x,t,Xtemp3,degree);
 }

 }

 c(g,f(y,x,t,g,degree))

} ######################################### akhir dari fungsi pso

#output = pso(N,itermaks,gcvmix,y,x,t,degree,nknot) ##################### Contoh cara panggil fungsi optimisasi PSO

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
	if (h == 0){Ktemp = 0}
	else {Ktemp = kernel(((x0-x[s])/h))}
	K = K + Ktemp
 }
 
 if (h == 0){Ktemp2 = 0}
 else {Ktemp2 = kernel(((x0-x1)/h))}
 
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


gcvmix = function(y,x,t,hknot,degree) #### hknot = [bandwidth | knot] (augmented)
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

#print(indeks[,1])

##### Data variabel x ########################

indeksx = removezero(indeks[,1]) ##### indeks variabel x ada di kolom ke-1 di file cvs indeks 
#print(length(indeksx))
x = matrix(ncol=length(indeksx),nrow=nrow(y))

for (i in 1:length(indeksx))
{
 for (j in 1:nrow(y))
 {
   x[j,i] = mydata[j,indeksx[i]+var_ind_start] #### variabel respon dimulai dari kolom ke-var_ind_start+1 di file csv data
 }
}
# print(x)
# print(ncol(x))


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

bandknotmgcv = pso(N,itermaks,gcvmix,y,x,t,batasvar,degree,nknot); ######## pencarian bandwidth optimum dengan PSO 

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
write.csv(output, "/Users/irwansyah/Documents/mixed-spline-project/experiment/result_a_pso_x1-5-kernel_x6-93-spline.csv")

parameters = c(degree,nknot,knotf,bandakhir);
write.csv(parameters, "/Users/irwansyah/Documents/mixed-spline-project/experiment/params_a_pso_x1-5-kernel_x6-93-spline.csv")

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
