# UNevada Reno Comparative Phylogenetics R Workshop
# Cheat Sheet
# Author: Hannah L. Owens, University of Copenhagen
# 30 July, 2019

# Libraries ----
library(ape); #The Swiss army knife of comparative phylogenetics R packages
library(geiger); #A package for "macroevolutionary simulation and estimating parameters related to diversification"
library(phytools); #Very useful for visualization particularly, great blog support
library(RColorBrewer); # Accessory package with better color
library(plotrix); #Contains a useful way of producing a color legend

# Loading Trees ----
istioTree <- read.nexus("Fish_12Tax_time_calibrated.tre"); # Can also load in newick files with read.tree(). 

## Looking at the tree
istioTree; #Shows the object
str(istioTree); #Shows how the object is structured

## How edges are structured
head(istioTree$edge); #Two columns, start node and end node/tip; tips labeled first
head(istioTree$edge.length); #Gives lengths of branches by edge

## How tips are structured, is it ultrametric
istioTree$tip.label; #Gives a list of tip labels, order corresponds to edge number
is.ultrametric(istioTree); #Checks to see if the tree is ultrametric

# Plotting Trees ----
plot(istioTree, cex = 2); #Basic view (phylogram)
nodelabels(cex = 2); #Shows numeric references for each node

# Searching Trees ----
## Search for taxa by name
colors <- c("black", "red");
matches <- grep("Kajikia", istioTree$tip.label, value = T); #Looking for taxa that contain "Kajikia"
plot(istioTree, cex = 2, 
     tip.col = ifelse(istioTree$tip.label %in% matches,'red','black')); #Plots taxa that contain "Kajikia" in red.

## Finding the node with a most-recent common ancestor.
matches; #Remember, these are our matches
findMRCA(tree = istioTree, tips = matches, type = "node"); #What is the node number of the most recent common ancestor?

## Some fancy tree tricks
plot(istioTree, type = "cladogram", use.edge.length = F, 
     edge.width = 2, font = 2); #Plotting a cladogram

plot(istioTree, type = "fan", rotate.tree = 30, 
     main = "Phylogeny of Billfishes", label.offset = 2,
     tip.color = c("blue", "red", "red", "yellow2", 
                   "yellow2", "orange1", "chartreuse4", 
                   "chartreuse4", "chartreuse4", 
                   "chartreuse4", "black", "purple"));#Plotting a fan with custom colored tips

## Extracting a clade of interest
subsetIstioTree <- extract.clade(istioTree, node = 18); # Tip: Nodes correspond to those plotted using nodelabels().


# Character Data ----
## Loading data
characterTable <- read.csv("CodingTableThresh95.csv", row.names = 1);
characterTable; #Showing the data

## Union with tree
treeWData <- treedata(istioTree, characterTable, sort = T); #Trims out taxa that are missing from the tree or the character table
plot(treeWData$phy, main = "Trimmed Tree", cex.main = 2, cex = 1); #Here's the trimmed tree
treeWData$data; #Here's the data

# Character evolution ----
#The basic function, in Geiger
charModel <- fitContinuous(treeWData$phy, treeWData$data[,17], model = "BM", SE = std.error(treeWData$data[,17]));
#How it works
charModel;
#How the results look
charModel$opt;
#The relevant stats

#Set up a table for results
reconAICcTable <- matrix(data = NA, nrow = 5, ncol = ncol(treeWData$data));
colnames(reconAICcTable) <- colnames(treeWData$data);
rownames(reconAICcTable) <- c("BM", "OU", "EB", "trend", "white");
reconRelWeightTable <- reconAICcTable; #Both tables will have the same structure
#Run through a loop that fits as set of character evolution models for each character
count <- 1;
while (count <= ncol(treeWData$data)){
  reconAICcTable[1,count] <- fitContinuous(treeWData$phy, treeWData$data[,count], model = "BM", 
                                           SE = std.error(treeWData$data[,count]))$opt$aicc; #Brownian motion
  reconAICcTable[2,count] <- fitContinuous(treeWData$phy, treeWData$data[,count], model = "OU", 
                                           SE = std.error(treeWData$data[,count]), bounds = list(alpha = c(min = exp(-500), max = exp(10))))$opt$aicc; #Ornstein-Uhlenbeck
  reconAICcTable[3,count] <- fitContinuous(treeWData$phy, treeWData$data[,count], model = "EB", 
                                           SE = std.error(treeWData$data[,count]), bounds = list(a = c(min = -3, max = 3)))$opt$aicc; #Early burst
  reconAICcTable[4,count] <- fitContinuous(treeWData$phy, treeWData$data[,count], model = "trend", 
                                           SE = std.error(treeWData$data[,count]))$opt$aicc; #Diffusion model w/ linear trend in rates through time
  reconAICcTable[5,count] <- fitContinuous(treeWData$phy, treeWData$data[,count], model = "white", 
                                           SE = std.error(treeWData$data[,count]))$opt$aicc; #White noise (non-phylogenetic model with no covariance among species)
  deltaAICc <- (reconAICcTable[,count])-min(reconAICcTable[,count]); #Calculates delta AICc for each model
  reconRelWeightTable[,count] <- exp(-.5*deltaAICc)/sum(exp(-.5*deltaAICc)); #Calculates relative weights of different models
  count <- count + 1;
}
round(reconRelWeightTable, digits = 4); #Displays relative weights of AICc, rounded to four digits

# Test phylogenetic signal----
#The basic functions
phytools::phylosig(treeWData$phy, x=treeWData$data[,1],method = "K");

#Set up a table for results
sigTable <- matrix(data = NA, ncol = 4, nrow = ncol(treeWData$data));
rownames(sigTable) <- colnames(treeWData$data);
colnames(sigTable) <- c("Lambda", "Lambda P-value", "K", "K P-value");
#Loop to calculate Blomberg's K and P-value for each character
count <- 1
while(count <= ncol(treeWData$data)){
  temp <- phytools::phylosig(treeWData$phy, treeWData$data[,count], 
                             method = "lambda", test = T); #Calculate Pagel's lambda
  sigTable[count,1] <- temp$lambda;
  sigTable[count,2] <- temp$P;
  temp <- phytools::phylosig(treeWData$phy, treeWData$data[,count], 
                             method = "K", test = T); #Calculate Blomberg's k
  sigTable[count,3] <- temp$K;
  sigTable[count,4] <- temp$P;
  count <- count + 1;
}
sigTable; #Strongest signal is for minimum annual SST, niche range

# Character Reconstructions: Discrete traits ----
discChar <- treeWData$data[,9]; #Getting your tip values

## Discrete trait reconstruction
recon <- ace(discChar, treeWData$phy, type = "discrete", 
               method = "ML", model = "ER"); #Do reconstruction.
recon$lik.anc; #What it looks like

## Plotting  
plot(treeWData$phy, main = "Reconstruction of preferred prey", 
     label.offset = 2.5, no.margin = F); #Plot tree
reconCol <- brewer.pal(3, "RdYlBu");   #Set up colors
nodelabels(text = rep(" ", 21-12+1), node = seq(12,21,1), 
           pie = recon$lik.anc, 
           piecol = c(reconCol[1], reconCol[2], reconCol[3])); #Labels nodes
tiplabels(text = rep(" ", 11), tip = seq(1,11,1), 
            frame = "circle", bg = reconCol[discChar]); #Labels tips
legend(x = 0, y = 5, fill = reconCol, bty = "n", 
       legend = c("Sardines", "Pompano", "Hamburgers")); #Plots legend

# Character Reconstructions: Continuous traits ----
## Do the reconstruction
recon <- ace(x = treeWData$data[,8], phy = treeWData$phy, 
             type = "continuous", method = "GLS", 
             corStruct = corBrownian(value=1, treeWData$phy));

## Create a matrix that contains each node's reconstructed values
reconMatrix <- matrix(
  nrow=(length(treeWData$phy$tip.label) + treeWData$phy$Nnode), 
  ncol=1);#Creating an empty matrix
colnames(reconMatrix) <- colnames(treeWData$data)[8]; #Labeling columns
rownames(reconMatrix) <- c(treeWData$phy$tip.label, 
                           seq(from = length(unlist(treeWData$phy$tip.label)) + 1, 
                               to = nrow(reconMatrix)));  #Labeling rows
reconMatrix[,1] <- round(c(treeWData$data[,8], recon$ace), 2); #Adding the tip data
reconMatrix <-as.matrix(reconMatrix); #Making sure it's in the right format

##Colors
colPal <- brewer.pal(9, "Reds"); #Definining the color palette
normScores <- round((reconMatrix[,1]-min(reconMatrix[,1]))/
                      (max(reconMatrix[,1])-min(reconMatrix[,1]))*8,0); #Normalizing the character scores on a scale of 1 to 9
reconCol <- colPal[normScores+1]; #Assigning colors to each character score

## Plotting
plot(treeWData$phy, main = "Reconstruction of niche breadth", 
     label.offset = 2.5, no.margin = F, cex = 1, cex.main = 2);
nodelabels(text = rep(" ", 21-12+1), node = seq(12,21,1), 
           cex = 0.65, frame = "circle", bg = reconCol[12:21]);
tiplabels(text = rep(" ", 11), seq(1,11,1), cex = 0.65, 
          frame = "circle", bg = reconCol[1:11]);

## Making a key for the plot
col.labels<-c("Narrow", "Broad")
color.legend(11,0,13,4,col.labels,
             colPal,gradient="y", cex = 1);

# Further Reading ----
## Analysis of Phylogenetics and Evolution with R
### By Emmanuel Paradis
#### Lead author of ape package
#### http://ape-package.ird.fr/APER.html
  
## Phytools Blog
### By Liam Revell
#### Lead author of Phytools package
#### http://blog.phytools.org/