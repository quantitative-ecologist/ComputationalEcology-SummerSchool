Tf <- 30 # game duration
########################  start the game in a patch   ########################
####################### WITHIN PATCH DECISION ################################

N_in_t0=50    # total number of marbles at a given station
N1_in_t0=10           #marbles of type 1
N2_in_t0=N1_in_t0                    #marbles of type2
Nsubstrat = N_in_t0-(N1_in_t0+N2_in_t0)
a_in = (60/10)                      # searching rate when not handling (marbles that can be drawn /min)
h1_in = (15/60)                      # handling time prey type 1 (min/marble)
h2_in = (15/60)                 # handling time prey type 1 (min/marble)
e1_in = 40           # points gained by consuming a prey of type 1 (points/marble) 
e2_in = 5                # points gained by consuming a prey of type 1 (points/marble)
p1_in=1                    #if marble of type 1 are taken, p = 1; 0 otherwise
p2_in=1              #if marble of type 2 are taken, p = 1; 0 otherwise

E_in_t0 <- ((a_in*(N1_in_t0/N_in_t0)*e1_in*p1_in)+(a_in*(N2_in_t0/N_in_t0)*e2_in*p2_in))/(1+(a_in*(N1_in_t0/N_in_t0)*h1_in*p1_in)+(a_in*(N2_in_t0/N_in_t0)*h2_in*p2_in))
cons1_in_T0 <- (a_in*(N1_in_t0/N_in_t0)*p1_in)/(1+(a_in*(N1_in_t0/N_in_t0)*h1_in*p1_in)+(a_in*(N2_in_t0/N_in_t0)*h2_in*p2_in))
cons2_in_T0 <- (a_in*(N2_in_t0/N_in_t0)*p2_in)/(1+(a_in*(N1_in_t0/N_in_t0)*h1_in*p1_in)+(a_in*(N2_in_t0/N_in_t0)*h2_in*p2_in))

E_in_t <- E_in_t0
#cumE_in_t <- E_in_t0
N1_in_t <- N1_in_t0
N2_in_t <- N2_in_t0
N_in_t <- N_in_t0
vN1_in_t <- numeric()
vN2_in_t <- numeric()
vN_in_t <- numeric()
vE_in_t <- numeric()
vnum <- numeric()
vdenom <- numeric()

########################################  
for(i in 1:(Tf)){
  N1_in_t <- N1_in_t-((a_in*(N1_in_t/N_in_t)*p1_in)/(1+(a_in*(N1_in_t/N_in_t)*h1_in*p1_in)+(a_in*(N2_in_t/N_in_t)*h2_in*p2_in)))
  N2_in_t <- N2_in_t-((a_in*(N2_in_t/N_in_t)*p2_in)/(1+(a_in*(N1_in_t/N_in_t)*h1_in*p1_in)+(a_in*(N2_in_t/N_in_t)*h2_in*p2_in)))
  E_in_t <- ((a_in*(N1_in_t/N_in_t)*e1_in*p1_in)+(a_in*(N2_in_t/N_in_t)*e2_in*p2_in))/(1+(a_in*(N1_in_t/N_in_t)*h1_in*p1_in)+(a_in*(N2_in_t/N_in_t)*h2_in*p2_in))
  N_in_t <- N1_in_t + N2_in_t + Nsubstrat
  num  <- ((a_in*(N1_in_t/N_in_t)*e1_in*p1_in)+(a_in*(N2_in_t/N_in_t)*e2_in*p2_in))
  denom <- (1+(a_in*(N1_in_t/N_in_t)*h1_in*p1_in)+(a_in*(N2_in_t/N_in_t)*h2_in*p2_in))
  vN1_in_t[i] <- (N1_in_t)
  vN2_in_t[i] <- (N2_in_t)
  vN_in_t[i] <- (N_in_t)
  vE_in_t[i] <- (E_in_t)
  vnum[i] <- (num)
  vdenom[i] <- (denom)
}
E_in <- vnum/vdenom  # rate of point collection 
E_in1 <- append((E_in_t0), E_in)
totgain_temp <- cumsum(E_in1)  # gain function (for each min)
totgain <- totgain_temp[1:Tf]  # gain function only during Tf  (to easily correct for poor coding!)

Residencytime <- 1:Tf
plot(Residencytime, totgain)

####################### AMONG PATCH DECISION ################################
a     <- 30*23                      #area search rate  (m2/min), i.e., width*speed 
N1    <- (200/(120000))  #  (200 patches; 30 participants, environ 12 ha; 120 000 m2, distance moy entre patches = 12.25 m (=1/(2*((200/120000)^0.5))) = 315 x 315 m)              ###### patch density  (patches/m2) if I put N1 = 1/a, N2 = 0.3/a, then 1 type is best instead of 2
N2    <- N1        ###### density of patch type 2  (patches/m2)
mus   <- 0.45                     ### survival probability when travelling between patches
muh   <- mus                   ### survival probability when at a patches (i.e., when handling food)
mur   <- .9                      ### survival probability when in refuge
p1 <- 1
p2 <- 1                 #p2=0 because only one type of resource patch

hx <- 1:Tf              # h: min spent in a patch handing marbles 
# encounter (patch/min) * handling a patch (min/patch)
Trx <- 1:Tf
hpatcha <- c(rep(hx, times=length(Trx)))
totgain1 <- c(rep(totgain, times=length(Trx)))
Tr <- c(rep(Trx, each=(length(hx))))

Alldata1 <- matrix(c(Tr, hpatcha, totgain1), nrow = (length(Tr)), ncol = 3)
colnames(Alldata1) <- c("Tr", "hpatch", "e_intra")

lambdaAll1 <- a*N1 #Encounter rate (m2/min*patch/m2 = patch/min) with patch of type 1 (only this type present)
# Minutes to encounter next patch = 1/lambdaAll1
Timepatch_enc1 <- 1/lambdaAll1 
timeto_search_handl1 <- Timepatch_enc1 + hpatcha # time spent in each patch + time to the next patch 
tot_patch_process_during_Tf <- (Tf-Tr)/timeto_search_handl1 #total no of patches visited during Tf
Num_patches<- tot_patch_process_during_Tf

remaining <- (Tf-Tr) %% timeto_search_handl1   # modulo to get min remaining 
diff1 <- remaining-hpatcha    # if +, then can consume the last patch completely, i.e., 2.7 patches --> 3 patches completely consumed
nb_patch_fullyhandl <- ifelse(diff1<0,floor(Num_patches), ceiling(Num_patches))  ### number of patch fully handled
bbbb <- ifelse(diff1<0,totgain_temp[remaining], 0)       #### gain in the last patch is not fully handled
cccc <- ifelse(diff1<0,0,diff1) ############ search time at the end... 0 if last patch not fully consumed, remaining otherwise
dddd <- ifelse(diff1<0,remaining,0)      #### time spent handling last patch if only partially consumed

gain_duringTf <- (nb_patch_fullyhandl*totgain1)+bbbb
Ts_inter <- (Timepatch_enc1*(floor(Num_patches)))+cccc     ### assume that consumer starts in a patch
Th1 <- (hpatcha*nb_patch_fullyhandl)+dddd

mu <- ((Ts_inter*mus+Th1*muh+Tr*mur)/Tf)   # overall survival probability, given the time spent in different behaviours (searching, handling marbles at a station, waiting in a refuge)
Eall1 <- gain_duringTf*mu; # Expected long term gain rate, given probability of survival

loc1 <- which(Eall1==max(Eall1))

zz <- matrix(c(Eall1,gain_duringTf, mu,hpatcha,nb_patch_fullyhandl, Ts_inter,Th1,Tr), ncol=8)
colnames(zz) <- c("Eall", "E", "mu", "h", "NbPatch_consumed", "Ts","Th","Tr")

cc1 <- c(Eall1[loc1],gain_duringTf[loc1], mu[loc1],hpatcha[loc1],nb_patch_fullyhandl[loc1], Ts_inter[loc1],Th1[loc1], Tr[loc1])
### Long-term expected gain, overall risk, handling rate, number of stations fully consumed, and time spent searching, handling, in refuge
cc1

timex <- 1:length(totgain_temp)
timex <- timex+Ts_inter[loc1]
gain <- totgain_temp/timex
plot(timex, gain)


