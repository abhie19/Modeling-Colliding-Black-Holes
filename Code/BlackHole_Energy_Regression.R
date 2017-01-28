# Import data files
data1 = read.table("/Users/Abhishek/Library/Mobile Documents/com~apple~CloudDocs/IU/SEM 3/Astronomy Data Analysis/Data/m10vx1.dat", header=FALSE)
data2 = read.table("/Users/Abhishek/Library/Mobile Documents/com~apple~CloudDocs/IU/SEM 3/Astronomy Data Analysis/Data/m10vx1p5.dat", header=FALSE)
data3 = read.table("/Users/Abhishek/Library/Mobile Documents/com~apple~CloudDocs/IU/SEM 3/Astronomy Data Analysis/Data/m10vx2.dat", header=FALSE)
data4 = read.table("/Users/Abhishek/Library/Mobile Documents/com~apple~CloudDocs/IU/SEM 3/Astronomy Data Analysis/Data/m10v75.dat", header=FALSE)
data5 = read.table("/Users/Abhishek/Library/Mobile Documents/com~apple~CloudDocs/IU/SEM 3/Astronomy Data Analysis/Data/New Data/m10vx1i2.dat", header=FALSE)
data6 = read.table("/Users/Abhishek/Library/Mobile Documents/com~apple~CloudDocs/IU/SEM 3/Astronomy Data Analysis/Data/New Data/m10vx1p5i2.dat", header=FALSE)
data7 = read.table("/Users/Abhishek/Library/Mobile Documents/com~apple~CloudDocs/IU/SEM 3/Astronomy Data Analysis/Data/New Data/m10vx2i2.dat", header=FALSE)
data8 = read.table("/Users/Abhishek/Library/Mobile Documents/com~apple~CloudDocs/IU/SEM 3/Astronomy Data Analysis/Data/New Data/m10vx1p75.dat", header=FALSE)
data9 = read.table("/Users/Abhishek/Library/Mobile Documents/com~apple~CloudDocs/IU/SEM 3/Astronomy Data Analysis/Data/New Data/m10vx1p25.dat", header=FALSE)

head(data9)

# Rename the headers
names(data1) <- c("Time", "Radius", "V3", "Energy", "V5","V6", "Orbit_Amplitude", "V8","Orbital_Period", "V10","V11","V12","V13","Velocity")
names(data2) <- c("Time", "Radius", "V3", "Energy", "V5","V6", "Orbit_Amplitude", "V8","Orbital_Period", "V10","V11","V12","V13","Velocity")
names(data3) <- c("Time", "Radius", "V3", "Energy", "V5","V6", "Orbit_Amplitude", "V8","Orbital_Period", "V10","V11","V12","V13","Velocity")
names(data4) <- c("Time", "Radius", "V3", "Energy", "V5","V6", "Orbit_Amplitude", "V8","Orbital_Period", "V10","V11","V12","V13","Velocity")
names(data5) <- c("Time", "Radius", "V3", "Energy", "V5","V6", "Orbit_Amplitude", "V8","Orbital_Period", "V10","V11","V12","V13","Velocity")
names(data6) <- c("Time", "Radius", "V3", "Energy", "V5","V6", "Orbit_Amplitude", "V8","Orbital_Period", "V10","V11","V12","V13","Velocity")
names(data7) <- c("Time", "Radius", "V3", "Energy", "V5","V6", "Orbit_Amplitude", "V8","Orbital_Period", "V10","V11","V12","V13","Velocity")
names(data8) <- c("Time", "Radius", "V3", "Energy", "V5","V6", "Orbit_Amplitude", "V8","Orbital_Period", "V10","V11","V12","V13","Velocity")
names(data9) <- c("Time", "Radius", "V3", "Energy", "V5","V6", "Orbit_Amplitude", "V8","Orbital_Period", "V10","V11","V12","V13","Velocity")

attach(data1)
attach(data2)
attach(data3)
attach(data4)
attach(data5)
attach(data6)
attach(data7)
attach(data8)
attach(data9)

# Explore Data

plot(x = data1_mod$Time, y=data1_mod$Energy, xlab = "Time", ylab="Energy", main = "Time vs Energy Data 1_mod") 
plot(x = data3_mod$Time, y=data3_mod$Energy, xlab = "Time", ylab="Radius", main = "Time vs Energy Data 9")
plot(x = data1$Time, y=data1$Orbit_Amplitude, xlab = "Time", ylab="Orbital Amplitude")
plot(x = data1$Time, y=data1$Orbital_Period, xlab = "Time", ylab="Orbital Period")

pairs(~Time+Energy+Radius+Orbit_Amplitude+Orbital_Period, data= data5, main ="Scatterplot Mix on Data 3")

# library(dplyr)
data1_mod <- subset(data1, data1$Time<=65)
data1_mod<-select(data1_mod, Time,Radius, Energy, Orbit_Amplitude, Orbital_Period, Velocity)

data2_mod <- subset(data2, data2$Time<=100)
data2_mod<-select(data2_mod, Time,Radius, Energy, Orbit_Amplitude, Orbital_Period)

data3_mod <- subset(data3, data3$Time<=225)
data3_mod<-select(data3_mod, Time,Radius, Energy, Orbit_Amplitude, Orbital_Period, Velocity)

data4_mod <- subset(data4, data4$Time<=300)
data4_mod<-select(data4_mod, Time,Radius, Energy, Orbit_Amplitude, Orbital_Period, Velocity)

data5_mod <- subset(data5, data5$Time<=150)
data5_mod<- select(data5_mod, Time,Radius, Energy, Orbit_Amplitude, Orbital_Period, Velocity)

data6_mod <- subset(data6, data6$Time<=150)
data6_mod<-select(data6_mod, Time,Radius, Energy, Orbit_Amplitude, Orbital_Period, Velocity)

data7_mod <- subset(data7, data7$Time<=215)
data7_mod<-select(data7_mod, Time,Radius, Energy, Orbit_Amplitude, Orbital_Period, Velocity)

data8_mod <- subset(data8, data8$Time<=200)
data8_mod<-select(data8_mod, Time,Radius, Energy, Orbit_Amplitude, Orbital_Period, Velocity)

data9_mod <- subset(data9, data9$Time<=100)
data9_mod<-select(data9_mod, Time,Radius, Energy, Orbit_Amplitude, Orbital_Period, Velocity)


# Adding Decay Rate Column

# For Data3_mod
for(n in 1:nrow(data9_mod)){
  data9_mod[n,7] <- (data9_mod[n+1,3]-data9_mod[n,3])/(data9_mod[n+1,1]-data9_mod[n,1])
} 
data9_mod[400,7] = mean(data9_mod[1:399,7])
# data3_mod <- data3_mod[1:900,]
colnames(data9_mod) <- c("Time", "Radius", "Energy", "Orbital_Amplitude", "Orbital_Period", "Velocity", "Decay")
# For Data1_mod

for(n in 1:nrow(data1_mod)){
  data1_mod[n,7] <- (data1_mod[n+1,3]-data1_mod[n,3])/(data1_mod[n+1,1]-data1_mod[n,1])
} 
data1_mod[260,7] = mean(data1_mod[1:259,7])

colnames(data1_mod) <- c("Time", "Radius", "Energy", "Orbital_Amplitude", "Orbital_Period", "Velocity", "Decay")

# Merging all data frames 
all_data <- merge(data1_mod,data2_mod, all = TRUE)
all_data <- merge(all_data,data3_mod, all = TRUE)
all_data <- merge(all_data,data9_mod, all = TRUE)

# Plotting with ggplot2
{
library(ggplot2)
p <- ggplot(all_data, aes(Time,Energy,color = V7)) + geom_point()
p
}
#Plotting time vs energy
plot(data1_mod$Time, -abs(data1_mod$Energy))

# Performing Log transformations on data
{
data1_mod["log_Energy"]<--(log(abs(data1_mod$Energy)))
}
# Creating linear model
{
fit1 <- lm(Energy ~ Time + Radius + Orbit_Amplitude + Orbital_Period, data=data1_mod)
fit2 <- lm(Energy ~ Time + Orbit_Amplitude + Orbital_Period, data=data1_mod)
 
summary(fit1)
abline(lm(fit1))

#Making test data
test_data <- select(data2_mod, Time, Radius, Orbit_Amplitude, Orbital_Period)
test_data <- test_data[1:20,]
}

# Code for taking only 1 standard deviation values
{
`%between%`<-function(x,rng) x>rng[1] & x<rng[2]
sd_data <- sd(data3_mod[,7])
mean_data <- mean(data3_mod[,7])
rng1 = mean_data-(1/2*sd_data)
rng2 = mean_data+(1/2*sd_data)
data3_mod_testsd <- data3_mod["Decay"] %between%  c(rng1,rng2)
data3_mod_filter <- data3_mod[data3_mod_testsd,]
plot(data3_mod_filter$Time,data3_mod_filter$Energy)
}
# Testing

library(ggplot2)
p <- ggplot(data3_mod_filter, aes(Time, Energy)) + geom_point()
# p <- ggplot(merge(data3_mod,data3_mod_filter, all= TRUE), aes(Time, Energy)) + geom_point()
p

# NLS PRACTICE START
{
# generate data
beta <- 0.05
n <- 100
temp <- data.frame(y = exp(beta * seq(n)) + rnorm(n), x = seq(n))

# plot data
plot(temp$x, temp$y)

# fit non-linear model
mod <- nls(y ~ exp(a + b * x), data = temp, start = list(a = 0, b = 0))

##################### WORKING 
m<-nls(y~a*exp(a*x), data = df, start= c(a = -0.004))

# add fitted curve
lines(temp$x, predict(mod, list(x = temp$x)))
}
# NLS PRACTICE END

# describing buckets and bucket values 
{
# 1st bucket Time 0 to 70 
new_energy = list[]
lower_range =0 
upper_range = 70
avg_decay1 = mean(data3_mod[data3_mod$Time>lower_range & data3_mod$Time<upper_range,7])
upper_index = nrow(data3_mod[data3_mod$Time>lower_range & data3_mod$Time<upper_range,])
lower_index = 1
# new_energy =list

initial_energy = data3_mod[lower_index,3] # set inital energy value to the first value
new_energy[lower_index] = initial_energy
for( i in (lower_index+1):upper_index){
  new_energy[i] = avg_decay1 * 0.25 + new_energy[i-1]
}
# Actual at 280 -0.6120911 ;Predicted at 280 -0.6121334
plot(data3_mod[lower_index:upper_index,1],new_energy[lower_index:upper_index])
# 2nd bucket Time 70 to 155
lower_index = upper_index
lower_range =70 
upper_range = 155
avg_decay2 = mean(data3_mod[data3_mod$Time>lower_range & data3_mod$Time<upper_range,7])
upper_index = nrow(data3_mod[data3_mod$Time>lower_range & data3_mod$Time<upper_range,])
upper_index = lower_index + upper_index


# initial_energy = data3_mod[lower_index,3] # set inital energy value to the first value
for( i in (lower_index+1):upper_index){
  new_energy[i] = avg_decay1 * 0.25 + new_energy[i-1]
}
# Actual at 620 -1.075842 ;Predicted at 620 -0.8655019 
plot(data3_mod[lower_index:upper_index,1],new_energy[lower_index:upper_index])
# avg_decay2 = mean(data3_mod[data3_mod$Time>70 & data3_mod$Time<155,7])
# 3rd bucket Time 155 to 250
lower_index = upper_index
lower_range =155 
upper_range = 250
avg_decay2 = mean(data3_mod[data3_mod$Time>lower_range & data3_mod$Time<upper_range,7])
upper_index = nrow(data3_mod[data3_mod$Time>lower_range & data3_mod$Time<upper_range,])
upper_index = lower_index + upper_index


# initial_energy = data3_mod[lower_index,3] # set inital energy value to the first value
# new_energy[lower_index] = initial_energy
for( i in (lower_index+1):upper_index){
  new_energy[i] = avg_decay1 * 0.25 + new_energy[i-1]
}
# Actual at 900 -1.450388 ;Predicted at 900 -1.074158 
plot(data3_mod[lower_index:upper_index,1],new_energy[lower_index:upper_index])
}

######### IMPORTANT FUNCTION ###############
# Automating process of orbit detection

# Create a new data frame to take care of orbit calculations FOR ANY DATASET
dataset = data9_mod
totalcount = dim(dataset)[1]
totalcount = totalcount-1

temp_df <- data.frame(matrix(c(1,2,3,4),ncol=4))
names(temp_df)<- c("Time", "Energy", "PreviousEnergy", "Decay")

i=0
direction = 0
prev_change ="TRUE"
prev_energy = dataset[1,3]
curr_energy = dataset[1,3]
j=1
for(i in (1:totalcount)){
  change = dataset[i,2]<dataset[i+1,2]
  if(prev_change!=change){
    direction = direction+1
  }
  if(direction == 2){
    # print(direction)
    # print(i)
    # print(data3_mod[i,1])
    # print(data3_mod[i,3])
    curr_energy = dataset[i,3]
    decay = prev_energy - (curr_energy) 
    direction = 0
    
    temp_df[j,1] <- dataset[i,1]
    temp_df[j,2] <- dataset[i,3]
    temp_df[j,3] <- prev_energy
    temp_df[j,4] <- decay
    
    print("One Orbit complete")
    print(paste0("The time is ", dataset[i,1]))
    print(paste0("The Previous Energy is ", prev_energy))
    print(paste0("The Energy is ", dataset[i,3]))
    print(paste0("The Decay is ", decay))
    j=j+1
  }
  prev_energy = curr_energy
  prev_change = change
}
# Change name of your DF
OrbitDF9<-temp_df 


#Some Rough
{# 
# # For data2_mod
# OrbitDF2 <- data.frame(matrix(c(1,2,3,4),ncol=4))
# names(OrbitDF2)<- c("Time", "Energy", "PreviousEnergy", "Decay")
# # Automating process of orbit detection
# i=0
# direction = 0
# prev_change ="TRUE"
# prev_energy = -0.9653524
# curr_energy = -0.9653524
# j=1
# for(i in (1:399)){
#   change = data2_mod[i,2]<data2_mod[i+1,2]
#   if(prev_change!=change){
#     direction = direction+1
#   }
#   if(direction == 2){
#     
#     curr_energy = data2_mod[i,3]
#     decay = prev_energy - (curr_energy) 
#     direction = 0
#     
#     OrbitDF2[j,1] <- data2_mod[i,1]
#     OrbitDF2[j,2] <- data2_mod[i,3]
#     OrbitDF2[j,3] <- prev_energy
#     OrbitDF2[j,4] <- decay
#     
#     j=j+1
#   }
#   prev_energy = curr_energy
#   prev_change = change
# }
# 
# # For data1_mod
# OrbitDF3 <- data.frame(matrix(c(1,2,3,4),ncol=4))
# names(OrbitDF3)<- c("Time", "Energy", "PreviousEnergy", "Decay")
# # Automating process of orbit detection
# i=0
# direction = 0
# prev_change ="TRUE"
# prev_energy = -1.366778
# curr_energy = -1.366778
# j=1
# for(i in (1:259)){
#   change = data1_mod[i,2]<data1_mod[i+1,2]
#   if(prev_change!=change){
#     direction = direction+1
#   }
#   if(direction == 2){
#     
#     curr_energy = data1_mod[i,3]
#     decay = prev_energy - (curr_energy) 
#     direction = 0
#     
#     OrbitDF3[j,1] <- data1_mod[i,1]
#     OrbitDF3[j,2] <- data1_mod[i,3]
#     OrbitDF3[j,3] <- prev_energy
#     OrbitDF3[j,4] <- decay
#     
#     j=j+1
#   }
#   prev_energy = curr_energy
#   prev_change = change
# }
}

# Combining all the Orbit DFs
OrbitDF_Combined_all <- rbind(as.matrix(OrbitDF),as.matrix(OrbitDF2),as.matrix(OrbitDF3),as.matrix(OrbitDF5),as.matrix(OrbitDF6),as.matrix(OrbitDF7),as.matrix(OrbitDF8),as.matrix(OrbitDF9))
OrbitDF_Combined_all <- as.data.frame(OrbitDF_Combined_all)
colnames(OrbitDF_Combined_all) <- c("Time","Energy","PreviousEnergy","Decay")
plot(OrbitDF_Combined_all$Energy,OrbitDF_Combined_all$Decay)

# Plottin Decay
library(ggplot2)
plot_combined <- ggplot( geom_point(data = newDF) + geom_point(data = OrbitDF_Combined) )
plot_combined

{# Draw a regression line on combined data
abline(lm(OrbitDF_Combined$Decay~OrbitDF_Combined$Energy))
# Taking only positive decays and drawing a regression line on combined data
OrbitDF_Combined_pos<-OrbitDF_Combined[OrbitDF_Combined$Decay>0,]
plot(OrbitDF_Combined_pos$Energy,OrbitDF_Combined_pos$Decay)
}

# Prediction

#modellm <- lm(Decay~PreviousEnergy, data= OrbitDF_Combined_all)
modellm <- lm(Decay ~ poly(PreviousEnergy,3), data= OrbitDF_Combined_all)
summary(modellm)
#abline((modellm))
predicted.intervals <- predict(modellm,data.frame(PreviousEnergy=OrbitDF_Combined_all$PreviousEnergy),interval='confidence',level=0.99)

plot(OrbitDF_Combined_all$PreviousEnergy,OrbitDF_Combined_all$Decay,col='deepskyblue4',xlab='Energy',ylab='Decay',main='Observed data')
lines(OrbitDF$Energy,predicted.intervals[1:68,1],col='red',lwd=1)

# MODEL WITH NLS
modelnls <- nls(Decay ~ a + exp(b * PreviousEnergy), data = OrbitDF_Combined_all, start = list(a = -0.9743, b = 1))
summary(modelnls)
plot(OrbitDF_Combined2$PreviousEnergy,OrbitDF_Combined2$Decay,col='deepskyblue4',xlab='q',main='Observed data')
lines(OrbitDF_Combined_all$PreviousEnergy,predicted.intervals[],col='red',lwd=1)

newDF = {}
decayDF ={}
new_energy = -0.9653524
# predict(modellm, data.frame(PreviousEnergy=c(-0.4042221)))
j=0
for(i in (1:55)){
  decay = predict(modellm, data.frame(PreviousEnergy=c(new_energy)))
  decayDF[j]<-decay
  new_energy<- new_energy-decay
  newDF[j] <- new_energy
  j = j+1
}
newDF[length(newDF)+1]<-newDF[length(newDF)]
plot(OrbitDF2$Time,newDF,col="red",xlab="Time", ylab="Predicted Energy")
points(OrbitDF2$Time,OrbitDF2$PreviousEnergy, col="blue")

## Pretty Plot for Decay vs Energy

plot(OrbitDF_Combined_all$PreviousEnergy,OrbitDF_Combined_all$Decay,col='deepskyblue4',xlab='Energy', ylab='Decay', main='Observed data')


names(newDF)<-c("Energy")
OrbitDF_Combined2<- OrbitDF_Combined
plot(OrbitDF_Combined2$Time,OrbitDF_Combined2$Energy)
plot(OrbitDF_Combined$Time,OrbitDF_Combined$Energy)
