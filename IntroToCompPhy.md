<style>
.reveal h1, .reveal h2, .reveal h3 {
  word-wrap: normal;
  -webkit-hyphens: none;
  -moz-hyphens: none;
  hyphens: none;
}

.section .reveal .state-background{
  background-color: gray;
}

</style>

Comparative Phylogenetics in R
========================================================
date: 28 January, 2020
author: Hannah L. Owens
autosize: true
width: 1400
height: 1000
font-family: 'Futura'

<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/0/0f/Anchovy_closeup.jpg/800px-Anchovy_closeup.jpg" style="width:50%; align:right;"></img>

Twitter Policy
========================================================
   
<img src="https://upload.wikimedia.org/wikipedia/en/thumb/9/9f/Twitter_bird_logo_2012.svg/1280px-Twitter_bird_logo_2012.svg.png" style="background-color:transparent; border:0px; box-shadow:none;"></img>
***
    
- Tweets are ok
- Attendance implies consent (wide-angle group shots)
- Asking permission before posting closer shots is just part of being a cool human.

Introduction
========================================================
type: section


Introduction: The Plan
========================================================
  
- **Introduction**    

Introduction: What R Cannot Do (On Its Own)
========================================================
incremental: true
- Phylogenetic inference
- Alignments
- BUT
- "ips" is a wrapper for all your favorites
  + raxml
  + mrbayes
  + beauti and BEAST
  + BayesTraits
  + TreeAnnotator
  + mafft and prank
  + gblocks and aliscore

***
<img src="https://freesvg.org/img/Hammer.png" style="width: 100%"></img>


Introduction: What R CAN Do!
========================================================
incremental: true
- Plotting
- Comparative phylogenetics
  + Modeling character evolution
  + Calculating phylogenetic signal
  + Ancestral character reconstruction
  + [And MUCH more](https://cran.r-project.org/web/views/Phylogenetics.html)!
    

***
<img src="https://media.giphy.com/media/DgLsbUL7SG3kI/giphy.gif" style="width: 100%;"></img>


Activity Dataset
========================================================
type: section

Activities: Anchovies
========================================================
   
<div align="left">
<img src="anchovy.jpg" width=100%>
</div>
***
   
- Bloom *et al*. 2018. Evolution of body size and trophic position in migratory fishes: a phylogenetic comparative analysis of Clupeiformes (anchovies, herring, shad, and allies). *Biol. J. Linn. Soc.* 125: 302 - 314.
- Body size and trophic level from FishBase
- Anadromy from previous work, updated w/ new info

Activities: Anchovies
========================================================
   
<div align="left">
<img src="anchovy.jpg" width=100%>
</div>
***
   
- Range from 2 to 40 cm
- Some freshwater, brackish, and marine
- Diadromous
  + Migrate between fresh and marine waters


The Plan
====================================
type: section

- Introduction    
- **Loading and Plotting Phylogenies**
- Character Data
- Modeling Character Evolution
- Calculating Phylogenetic Signal
- Ancestral Character State Reconstruction

Loading and Plotting Phylogenies
====================================
type: section

Loading and Plotting Phylogenies
========================================================

- Loading tree files
- Parts of tree files
- Displaying trees



Loading and Plotting Phylogenies
========================================================

The libraries we're using


```r
library(ape)
  #The Swiss army knife of comparative phylogenetics R packages
library(geiger)
  #A package for "macroevolutionary simulation and estimating parameters related to diversification"
library(phytools)
  #Very useful for visualization particularly, great blog support
```

Loading and Plotting Phylogenies : Loading Tree
========================================================

- Can read in a single tree or multiple trees.
- Takes ".tre" OR ".nex" nexus files
- You can also load in newick files with read.tree(). 


```r
istioTree <- read.nexus("Data/Fish_12Tax_time_calibrated.tre")
```

Loading and Plotting Phylogenies : Parts of tree files
========================================================


```r
istioTree
```

```

Phylogenetic tree with 12 tips and 11 internal nodes.

Tip labels:
	Istiompax_indica, Istiophorus_albicans, Istiophorus_platypterus, Kajikia_albida, Kajikia_audax, Makaira_nigricans, ...

Rooted; includes branch lengths.
```

```r
str(istioTree)
```

```
List of 4
 $ edge       : int [1:22, 1:2] 13 14 15 16 17 17 16 15 18 19 ...
 $ edge.length: num [1:22] 19.98 39.23 6.64 5.47 8.79 ...
 $ Nnode      : int 11
 $ tip.label  : chr [1:12] "Istiompax_indica" "Istiophorus_albicans" "Istiophorus_platypterus" "Kajikia_albida" ...
 - attr(*, "class")= chr "phylo"
 - attr(*, "order")= chr "cladewise"
```

Loading and Plotting Phylogenies : Parts of tree files
========================================================

What does the file look like?


```r
head(istioTree$edge) #Two columns, start node and end node/tip; tips labeled first
```

```
     [,1] [,2]
[1,]   13   14
[2,]   14   15
[3,]   15   16
[4,]   16   17
[5,]   17    3
[6,]   17    2
```

```r
head(istioTree$edge.length) #Gives lengths of branches by edge
```

```
[1] 19.976147 39.227748  6.635844  5.473333  8.786535  8.786535
```

Loading and Plotting Phylogenies : Parts of tree files
========================================================

What does the file look like?


```r
istioTree$tip.label #Gives a list of tip labels, order corresponds to edge number
```

```
 [1] "Istiompax_indica"           "Istiophorus_albicans"      
 [3] "Istiophorus_platypterus"    "Kajikia_albida"            
 [5] "Kajikia_audax"              "Makaira_nigricans"         
 [7] "Tetrapturus_angustirostris" "Tetrapturus_belone"        
 [9] "Tetrapturus_georgii"        "Tetrapturus_pfluegeri"     
[11] "Trachurus_trachurus"        "Xiphias_gladius"           
```

```r
is.ultrametric(istioTree) #Checks to see if the tree is ultrametric
```

```
[1] TRUE
```

Loading and Plotting Phylogenies : Displaying trees
========================================================

The basics


```r
plot(istioTree, cex = 4) #Basic view (phylogram)
nodelabels(cex = 4)
```

<img src="IntroToCompPhy-figure/unnamed-chunk-6-1.png" title="plot of chunk unnamed-chunk-6" alt="plot of chunk unnamed-chunk-6" style="display: block; margin: auto;" />

Loading and Plotting Phylogenies : Searching trees
========================================================

Finding tips that match  
- not specific to comparative phylogenetics  


```r
colors <- c("black", "red")
matches <- grep("Kajikia", istioTree$tip.label, value = T)
matches
```

```
[1] "Kajikia_albida" "Kajikia_audax" 
```

Loading and Plotting Phylogenies : Searching trees
========================================================

Highlighting matching taxa


```r
plot(istioTree, cex = 4, 
     tip.col = ifelse(istioTree$tip.label %in% matches,
                      'red','black'))
```

<img src="IntroToCompPhy-figure/unnamed-chunk-8-1.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" style="display: block; margin: auto;" />

Loading and Plotting Phylogenies : Searching trees
========================================================

Finding the node with a most-recent common ancestor.


```r
#Remember, these are our matches
matches
```

```
[1] "Kajikia_albida" "Kajikia_audax" 
```

```r
#What is the node number?
findMRCA(tree = istioTree, tips = matches, type = "node")
```

```
[1] 20
```

Challenge 1: Load and Display a Tree
========================================================
type: prompt

We have given you a phylogeny of clupeiform fishes from 
[Bloom and Lovejoy, 2014](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2013.2081). Load it into R and display the tree. The file is called "Clupeiform phylogeny.nex".

Bonus challenge: Which node number represents the most
recent common ancestor for the *Anchoa* clade?

Challenge 1: Load and Display a Tree
========================================================
type: prompt

The node of the *Anchoa* MRCA is 178.


```r
anchovyTree <- read.nexus("Data/Clupeiform phylogeny.nex")
matches <- grep("Anchoa", anchovyTree$tip.label, value = T)
plot(anchovyTree, cex = 0.5, 
     tip.col = ifelse(anchovyTree$tip.label %in% matches,'red','black'))
findMRCA(tree = anchovyTree, tips = matches, type = "node")
```

```
[1] 178
```

```r
nodelabels(node = 178, text = "MRCA", bg = F, col = "red")
```

<img src="IntroToCompPhy-figure/unnamed-chunk-10-1.png" title="plot of chunk unnamed-chunk-10" alt="plot of chunk unnamed-chunk-10" style="display: block; margin: auto;" />

Loading and Plotting Phylogenies : Displaying trees
========================================================

Getting fancy...


```r
plot(istioTree, type = "cladogram", use.edge.length = F, 
     edge.width = 2, font = 4, cex = 2.5)
```

<img src="IntroToCompPhy-figure/unnamed-chunk-11-1.png" title="plot of chunk unnamed-chunk-11" alt="plot of chunk unnamed-chunk-11" style="display: block; margin: auto;" />

Loading and Plotting Phylogenies : Displaying trees
========================================================

Getting fancy...


```r
plot(istioTree, type = "fan", rotate.tree = 30, 
     main = "Phylogeny of Billfishes", cex.main = 3.5, 
     label.offset = 2, cex = 2.5, 
     tip.color = c("blue", "red", "red", "yellow2",
                   "yellow2", "orange1", "chartreuse4", 
                   "chartreuse4", "chartreuse4", 
                   "chartreuse4", "black", "purple"))
```

<img src="IntroToCompPhy-figure/unnamed-chunk-12-1.png" title="plot of chunk unnamed-chunk-12" alt="plot of chunk unnamed-chunk-12" style="display: block; margin: auto;" />

Challenge 2: Make Your Own "Fancy" Tree
========================================================
type: prompt

Try to imitate this tree. It is a subset of the anchovy tree trimmed using extract.clade() like this:
`subsetAnchovyTree <- extract.clade(anchovyTree, node = 160)`.
Use `?ape::plot.phylo()` for help.

<img src="IntroToCompPhy-figure/unnamed-chunk-13-1.png" title="plot of chunk unnamed-chunk-13" alt="plot of chunk unnamed-chunk-13" style="display: block; margin: auto;" />

Challenge 2: Make Your Own "Fancy" Tree
========================================================
type: prompt

Here's how I made that tree. How'd you do?


```r
subsetAnchovyTree <- extract.clade(anchovyTree, node = 160)
plot(subsetAnchovyTree, type = "fan", rotate.tree = 45, open.angle = 180, 
     main = "Phylogeny of Anchovies", cex.main = 5, col.main = "green",
     label.offset = 2, cex = 2.5, tip.color = c("blue", "red"))
```

<img src="IntroToCompPhy-figure/unnamed-chunk-14-1.png" title="plot of chunk unnamed-chunk-14" alt="plot of chunk unnamed-chunk-14" style="display: block; margin: auto;" />

The Plan
====================================
type: section

- Introduction    
- Loading and Plotting Phylogenies
- **Character Data**
- Modeling Character Evolution
- Calculating Phylogenetic Signal
- Ancestral Character State Reconstruction


Character Data
====================================
type: section

Character Data
====================================

- Loading data
- Uniting with tree
- Plotting on tree

Character Data : Loading data
========================================================

Multiple characters in a .csv file


```r
characterTable <- read.csv("Data/CodingTableThresh95.csv", row.names = 1)
characterTable[,1:4] #Showing the data
```

```
                           bathmin bathmax SSSmin SSS8max
Istiompax_indica             -7000    -485   3416    3651
Istiophorus_albicans         -6275      -1   3393    3741
Istiophorus_platypterus      -7272      -1   3359    3639
Kajikia_albida               -5750   -1267   3561    3681
Kajikia_audax                -7000      -1   3451    3638
Makaira_nigricans            -5750   -1107   3437    3694
Tetrapturus_angustirostris   -7000    -841   3425    3649
Tetrapturus_belone           -1950     -95   3677    3800
Tetrapturus_georgii          -6414   -3392   3458    3741
Tetrapturus_pfluegeri        -6693   -2112   3437    3741
Xiphias_gladius              -6738    -748   3337    3748
```

Character Data : Union with tree
========================================================

Combining the tree and the character data


```r
treeWData <- treedata(istioTree, characterTable, sort = T)
```

```
Warning in treedata(istioTree, characterTable, sort = T): The following tips were not found in 'data' and were dropped from 'phy':
	Trachurus_trachurus
```

```r
  #Trims out taxa that are missing from the tree or the character table
```

Character Data : Union with tree
========================================================

Original tree versus trimmed tree

![plot of chunk unnamed-chunk-17](IntroToCompPhy-figure/unnamed-chunk-17-1.png)![plot of chunk unnamed-chunk-17](IntroToCompPhy-figure/unnamed-chunk-17-2.png)

Character Data : Union with tree
========================================================

Combining the tree and the character data


```r
treeWData$data[,1:4] #Here's the data
```

```
                           bathmin bathmax SSSmin SSS8max
Istiompax_indica             -7000    -485   3416    3651
Istiophorus_albicans         -6275      -1   3393    3741
Istiophorus_platypterus      -7272      -1   3359    3639
Kajikia_albida               -5750   -1267   3561    3681
Kajikia_audax                -7000      -1   3451    3638
Makaira_nigricans            -5750   -1107   3437    3694
Tetrapturus_angustirostris   -7000    -841   3425    3649
Tetrapturus_belone           -1950     -95   3677    3800
Tetrapturus_georgii          -6414   -3392   3458    3741
Tetrapturus_pfluegeri        -6693   -2112   3437    3741
Xiphias_gladius              -6738    -748   3337    3748
```

Challenge 3: Uniting anchovy characters with the tree
========================================================
type: prompt

You've already read in the anchovy phylogeny. Now read in 
the character table, Bloom_etal_2018_Reduced_Dataset.csv,
from [Bloom *et al*., 2018](https://doi.org/10.1093/biolinnean/bly106).

How many taxa are represented in the phylogeny *and* the
character dataset?

Challenge 3: Uniting anchovy characters with the tree
========================================================
type: prompt


```r
anchovyCharacters <- read.csv(file = "Data/Bloom_etal_2018_Reduced_Dataset.csv", 
                              row.names = 1)
anchovyTreeWData <- treedata(anchovyTree, anchovyCharacters, sort = T);
plot(anchovyTreeWData$phy, main = "Anchovies!")
paste0("There are ", nrow(anchovyTreeWData$data), " taxa represented in both datasets.")
```

```
[1] "There are 50 taxa represented in both datasets."
```

![plot of chunk unnamed-chunk-19](IntroToCompPhy-figure/unnamed-chunk-19-1.png)

Plotting Characters on Trees
===============================================================
type: section

Plotting Characters on Trees
====================================
- Discrete traits
- Continuous traits

Plotting Characters on Trees: Discrete traits
===============================================================
Plotting discrete traits on tips.


```r
library(RColorBrewer) # Accessory package with better color

#Entering in discrete tip values
discChar <- c(1, 1, 1, 3, 3, 1, 2, 2, 2, 2, 1)

#Plotting
plot(treeWData$phy, main = "Preferred prey", cex.main = 4, 
     label.offset = 2.5, no.margin = F, cex = 3)

#Setting up colors.
reconCol <- brewer.pal(3, "RdYlBu")

#Plot tips
tiplabels(text = rep(" ", 11), tip = seq(1,11,1), cex = 3, 
          frame = "circle", bg = reconCol[discChar])

#Plot legend
legend(x = 0, y = 5, fill = reconCol, bty = "n", cex = 3, 
       legend = c("Sardines", "Pompano", "Hamburgers"))
```

Plotting Characters on Trees: Discrete traits
===============================================================
Plotting discrete traits on tips.

![plot of chunk unnamed-chunk-21](IntroToCompPhy-figure/unnamed-chunk-21-1.png)

Challenge 4: Plotting anchovy characters on the tree
========================================================
type: prompt

You have a tree united with characters. It's plotting time!

Your challenge: Color-code the tip labels according to whether a species is diadromous or not.


Challenge 4: Plotting anchovy characters on the tree
========================================================
type: prompt


```r
# Set up colors
labelCols <- c("red", "blue");
diadLabs <- labelCols[as.factor(anchovyTreeWData$data[,3])]

#Plotting
plot(anchovyTreeWData$phy, main = "Diadromy", cex.main = 4, 
     no.margin = F, cex = 3, 
     label.offset = 2.5)

tiplabels(text = rep(" ", length(anchovyTreeWData$phy$tip.label)), 
          cex = 1, frame = "circle", bg = diadLabs)

#Plot legend
legend(x = 0, y = 50, bty = "n", cex = 2, 
       fill = labelCols, 
       legend = c("Diadromous", "Non-diadromous")
```

Challenge 4: Plotting anchovy characters on the tree
========================================================
type: prompt

![plot of chunk unnamed-chunk-23](IntroToCompPhy-figure/unnamed-chunk-23-1.png)

Plotting Characters on Trees: Continuous traits
===============================================================
Plotting continuous characters using a color gradient


```r
library(plotrix); #Contains a useful way of producing a color legend

#Set up colors and assigning values
colPal <- brewer.pal(9, "Reds")
normedTips <- round((treeWData$data[,8]-min(treeWData$data[,8]))/
                           (max(treeWData$data[,8])-min(treeWData$data[,8]))*8,0)
reconCol <- colPal[normedTips+1]
reconCol
```

```
 [1] "#FB6A4A" "#FC9272" "#FC9272" "#FEE0D2" "#FC9272" "#FEE0D2" "#FB6A4A"
 [8] "#FFF5F0" "#FEE0D2" "#FCBBA1" "#67000D"
```


Plotting Characters on Trees: Continuous traits
===============================================================
Plotting continuous characters using a color gradient


```r
# Do the plot
plot(treeWData$phy, main = "Niche breadth", label.offset = 2.5, 
     no.margin = F)

#Label Tips
tiplabels(text = rep(" ", 11), seq(1,11,1), 
          cex = 0.65, frame = "circle", bg = reconCol[1:11])

#Making a key for the plot
col.labels<-c("Narrow", "Broad")
color.legend(11,0,13,4, col.labels, 
             colPal,gradient="y") #Numbers are plot coordinates
```

Plotting Characters on Trees: Continuous traits
===============================================================
Plotting continuous characters using a color gradient

![plot of chunk unnamed-chunk-26](IntroToCompPhy-figure/unnamed-chunk-26-1.png)

Challenge 5: Plotting continuous anchovy characters on the tree
========================================================
type: prompt

You guessed it. Time for a continuous character plot.

Your challenge: Color-code the tip labels according to log body size.

Bonus challenge: Add another tip label showing trophic position.
(Hint: check out `?tiplabels()`)

Challenge 5: Plotting continuous anchovy characters on the tree
========================================================
type: prompt


```r
# Set up colors
logBody <- as.numeric(anchovyTreeWData$data[,1])
colPalBody <- rev(brewer.pal(9, "Reds"))
normedTipsB <- round((logBody-min(logBody))/
                           (max(logBody)-min(logBody))*8,0)
tipColBody <- colPalBody[normedTipsB+1]

trophPos <- as.numeric(anchovyTreeWData$data[,2])
colPalTroph <- rev(brewer.pal(9, "Blues"))
normedTipsT <- round((trophPos-min(trophPos))/
                           (max(trophPos)-min(trophPos))*8,0)
tipColTroph <- colPalTroph[normedTipsT+1]
```

Challenge 5: Plotting continuous anchovy characters on the tree
========================================================
type: prompt


```r
#Plotting
plot(anchovyTreeWData$phy, main = "Body Size and Trophic level", 
     cex.main = 4, no.margin = F, cex = 3, label.offset = 5)

tiplabels(cex = 3, adj = 1, pch = 21, bg = tipColBody)

tiplabels(cex = 3, adj = 3.5, pch = 21, bg = tipColTroph)

#Making keys for the plot
col.labels<- round(c(min(
  as.numeric(anchovyTreeWData$data[,1])), 
  max(as.numeric(anchovyTreeWData$data[,1]))),2)
text(5,46, "Body Size", cex = 3)
color.legend(0,40,4,44,col.labels,rev(colPalBody),
             gradient="y", cex = 3, align = "rb")

col.labels<-round(c(min(
  as.numeric(anchovyTreeWData$data[,2])), 
  max(as.numeric(anchovyTreeWData$data[,2]))),2)
text(5,36, "Trophic Level", cex = 3)
color.legend(0,30,4,34,col.labels,rev(colPalTroph),
             gradient="y", cex = 3, align = "rb")
```

Challenge 5: Plotting continuous anchovy characters on the tree
========================================================
type: prompt

![plot of chunk unnamed-chunk-29](IntroToCompPhy-figure/unnamed-chunk-29-1.png)

The Plan
====================================
type: section

- Introduction    
- Loading and Plotting Phylogenies
- Character Data
- **Modeling Character Evolution**
- Calculating Phylogenetic Signal
- Ancestral Character State Reconstruction


Modeling Character Evolution
====================================
type: section

Modeling Character Evolution
====================================

- The basic function
- Setting up for the loop
- The loop
- The results

Modeling character evolution: Basic function
====================================


```r
library(plotrix) #Accessories for phylogenies
charModel <- fitContinuous(treeWData$phy, 
                           treeWData$data[,8], 
                           model = "BM", 
                           SE = std.error(treeWData$data[,8]));
charModel #How the results look
```

```
GEIGER-fitted comparative model of continuous data
 fitted 'BM' model parameters:
	sigsq = 43607453068390600.000000
	z0 = 2554614105.734349

 model summary:
	log-likelihood = -244.522632
	AIC = 493.045264
	AICc = 494.545264
	free parameters = 2

Convergence diagnostics:
	optimization iterations = 100
	failed iterations = 0
	frequency of best fit = 1.00

 object summary:
	'lik' -- likelihood function
	'bnd' -- bounds for likelihood search
	'res' -- optimization iteration summary
	'opt' -- maximum likelihood parameter estimates
```

Modeling character evolution: Basic function
====================================

The relevant stats from fitContinuous()


```r
charModel$opt
```

```
$sigsq
[1] 4.360745e+16

$z0
[1] 2554614106

$lnL
[1] -244.5226

$method
[1] "Brent"

$k
[1] 2

$aic
[1] 493.0453

$aicc
[1] 494.5453
```

Challenge 6: Anchovy character evolution
========================================================
type: prompt

Run fitContinuous() using two models, Brownian motion and early burst. Which model fits the data better (smallest AICc score).

NOTE: You may run into problems with the $data object if R converted your numbers to characters. Use "as.numeric()" and "names" to fix this!

Challenge 6: Anchovy character evolution
========================================================
type: prompt

Oh hey, look at that. Brownian motion is a better fit.


```r
# Data prep
bodySize <- as.numeric(anchovyTreeWData$data[,1])
names(bodySize) <- rownames(anchovyTreeWData$data)

# Fit models
chovyBMcharModel <- fitContinuous(anchovyTreeWData$phy, 
                           bodySize, 
                           model = "BM", 
                           SE = std.error(bodySize))
chovyEBModel <- fitContinuous(anchovyTreeWData$phy, 
                           bodySize, 
                           model = "EB", 
                           SE = std.error(bodySize))

paste0("AICc for BM: ", chovyBMcharModel$opt$aicc, ". AICc for EB: ", chovyEBModel$opt$aicc, ".")
```

```
[1] "AICc for BM: 17.1792349761662. AICc for EB: 19.4459022517872."
```

Modeling character evolution: Set up results table
====================================

- BM: Brownian motion
- OU: Ornstein-Uhlenbeck
- EB: Early burst
- trend: Diffusion model w/ linear trend in rates through time
- white: White noise (non-phylogenetic model with no covariance among species)


```r
reconAICcTable <- matrix(data = NA, nrow = 5, ncol = ncol(treeWData$data))
colnames(reconAICcTable) <- colnames(treeWData$data)
rownames(reconAICcTable) <- c("BM", "OU", "EB", "trend", "white")
reconRelWeightTable <- reconAICcTable; #Both tables have the same structure
```

Modeling character evolution: The loop
====================================

This loop fits a set of character evolution models for each character.


```r
count <- 1;
while (count <= ncol(treeWData$data)){
  reconAICcTable[1,count] <- fitContinuous(treeWData$phy, 
      treeWData$data[,count], model = "BM", 
      SE = std.error(treeWData$data[,count]))$opt$aicc
  reconAICcTable[2,count] <- fitContinuous(treeWData$phy, 
      treeWData$data[,count], model = "OU", 
      SE = std.error(treeWData$data[,count]), 
      bounds = list(alpha = c(min = exp(-500), 
      max = exp(10))))$opt$aicc;
  reconAICcTable[3,count] <- fitContinuous(treeWData$phy, 
      treeWData$data[,count], model = "EB", 
      SE = std.error(treeWData$data[,count]), 
      bounds = list(a = c(min = -3, max = 3)))$opt$aicc
  reconAICcTable[4,count] <- fitContinuous(treeWData$phy, 
      treeWData$data[,count], model = "trend", 
      SE = std.error(treeWData$data[,count]))$opt$aicc
  reconAICcTable[5,count] <- fitContinuous(treeWData$phy, 
      treeWData$data[,count], model = "white", 
      SE = std.error(treeWData$data[,count]))$opt$aicc
  deltaAICc <- (reconAICcTable[,count])-min(reconAICcTable[,count])
  reconRelWeightTable[,count] <- exp(-.5*deltaAICc)/
                                  sum(exp(-.5*deltaAICc))
  count <- count + 1
}
```

Modeling character evolution: The results
====================================

Relative weights of AICc, rounded to four digits:


```r
round(reconRelWeightTable, digits = 5)
```

```
      bathmin bathmax  SSSmin SSS8max SST13min SST13max tempRange nicheRange
BM    0.00094 0.00988 0.00672 0.02128  0.20715  0.00089   0.07440    0.45167
OU    0.10939 0.10807 0.10852 0.12026  0.08281  0.10939   0.10508    0.07090
EB    0.10938 0.10805 0.10852 0.12026  0.08281  0.10939   0.10508    0.07090
trend 0.00038 0.00350 0.00248 0.00709  0.04512  0.00037   0.02010    0.06776
white 0.77991 0.77050 0.77375 0.73112  0.58210  0.77996   0.69534    0.33876
         Diet
BM    0.64061
OU    0.11528
EB    0.11528
trend 0.12883
white 0.00000
```

Challenge 7: Anchovy character evolution
========================================================
type: prompt

What is the best model of character evolution for log body size and trophic position in anchovies?

Remember R changed all your numbers to characters. Grr.

Challenge 7: Anchovy character evolution
========================================================
type: prompt

Step 1: Set up the data and the table.


```r
# Data
chovyCharTab <- cbind(as.numeric(anchovyTreeWData$data[,1]),as.numeric(anchovyTreeWData$data[,2]))
rownames(chovyCharTab) <- rownames(anchovyTreeWData$data)
colnames(chovyCharTab) <- colnames(anchovyTreeWData$data)[1:2]

# Results tables
reconAICcTable <- matrix(data = NA, nrow = 5, ncol = ncol(chovyCharTab))
colnames(reconAICcTable) <- colnames(chovyCharTab)
rownames(reconAICcTable) <- c("BM", "OU", "EB", "trend", "white")
```

Challenge 7: Anchovy character evolution
========================================================
type: prompt

Step 2: Run that loop!


```r
count <- 1;
while (count <= ncol(chovyCharTab)){
  reconAICcTable[1,count] <- fitContinuous(anchovyTreeWData$phy, 
      chovyCharTab[,count], model = "BM", 
      SE = std.error(chovyCharTab[,count]))$opt$aicc;
  reconAICcTable[2,count] <- fitContinuous(anchovyTreeWData$phy, 
      chovyCharTab[,count], model = "OU", 
      SE = std.error(chovyCharTab[,count]), 
      bounds = list(alpha = c(min = exp(-500), 
      max = exp(10))))$opt$aicc;
  reconAICcTable[3,count] <- fitContinuous(anchovyTreeWData$phy, 
      chovyCharTab[,count], model = "EB", 
      SE = std.error(chovyCharTab[,count]), 
      bounds = list(a = c(min = -3, max = 3)))$opt$aicc
  reconAICcTable[4,count] <- fitContinuous(anchovyTreeWData$phy, 
      chovyCharTab[,count], model = "trend", 
      SE = std.error(chovyCharTab[,count]))$opt$aicc
  reconAICcTable[5,count] <- fitContinuous(anchovyTreeWData$phy, 
      chovyCharTab[,count], model = "white", 
      SE = std.error(chovyCharTab[,count]))$opt$aicc
  count <- count + 1
}
```

Challenge 7: Anchovy character evolution
========================================================
type: prompt

Step 3: Which model?


```r
reconAICcTable
```

```
      logbodysize trophic_position
BM       17.17923        -100.7962
OU       17.31522        -156.8098
EB       17.31522        -156.6843
trend    16.98079        -107.3646
white    49.82883        -159.0762
```

The Plan
====================================
type: section

- Introduction    
- Loading and Plotting Phylogenies
- Character Data
- Modeling Character Evolution
- **Calculating Phylogenetic Signal**
- Ancestral Character State Reconstruction


Calculating Phylogenetic Signal
====================================
type: section

Calculating Phylogenetic Signal
====================================

- The basic function
- Setting up for the loop
- The loop
- The results

Phylogenetic signal: The basic function
===============================================================


```r
phytools::phylosig(treeWData$phy, x=treeWData$data[,1],method = "K")
```

```
[1] 0.1342713
```

Phylogenetic signal: Setting up for the loop
===============================================================


```r
sigTable <- matrix(data = NA, ncol = 4, nrow = ncol(treeWData$data))
rownames(sigTable) <- colnames(treeWData$data)
colnames(sigTable) <- c("Lambda", "Lambda P-value", "K", "K P-value")
```

Phylogenetic signal: The loop
===============================================================
Pagel's Lambda, Blomberg's K, and their respective P-values for each character


```r
count <- 1
while(count <= ncol(treeWData$data)){
  temp <- phytools::phylosig(treeWData$phy, treeWData$data[,count], 
                             method = "lambda", test = T); #Calculate Pagel's lambda
  sigTable[count,1] <- temp$lambda
  sigTable[count,2] <- temp$P
  temp <- phytools::phylosig(treeWData$phy, treeWData$data[,count], 
                             method = "K", test = T); #Calculate Blomberg's k
  sigTable[count,3] <- temp$K
  sigTable[count,4] <- temp$P
  count <- count + 1
}
```

Phylogenetic signal: The results
===============================================================
Et voila!


```r
sigTable
```

```
                 Lambda Lambda P-value         K K P-value
bathmin    6.901244e-05   1.000000e+00 0.1342713     0.835
bathmax    6.901244e-05   1.000000e+00 0.1958716     0.700
SSSmin     6.754471e-05   1.000000e+00 0.2386497     0.558
SSS8max    6.901244e-05   1.000000e+00 0.2418311     0.546
SST13min   6.574920e-01   5.737320e-01 0.4547630     0.100
SST13max   6.901244e-05   1.000000e+00 0.1375421     0.820
tempRange  5.217411e-05   1.000000e+00 0.4011244     0.163
nicheRange 7.364686e-01   1.068847e-01 0.8079658     0.060
Diet       1.043882e+00   7.949707e-05 1.1791841     0.002
```

Challenge 8: Anchovy character evolution
========================================================
type: prompt

What is the phylogenetic signal for log body size and trophic position in anchovies?

Challenge 8: Anchovy phylogenetic signal
========================================================
type: prompt

Step 1: Set up table for results. 

```r
#Set up a table for results
sigTable <- matrix(data = NA, ncol = 4, nrow = ncol(chovyCharTab))
rownames(sigTable) <- colnames(chovyCharTab)
colnames(sigTable) <- c("Lambda", "Lambda P-value", "K", "K P-value")
```

Challenge 8: Anchovy phylogenetic signal
========================================================
type: prompt

Step 2: Loop through each column. 

```r
#Loop to calculate Blomberg's K and P-value for each character
count <- 1
while(count <= ncol(chovyCharTab)){
  temp <- phytools::phylosig(anchovyTreeWData$phy, chovyCharTab[,count], 
                             method = "lambda", test = T); #Calculate Pagel's lambda
  sigTable[count,1] <- temp$lambda
  sigTable[count,2] <- temp$P
  temp <- phytools::phylosig(anchovyTreeWData$phy, chovyCharTab[,count], 
                             method = "K", test = T) #Calculate Blomberg's k
  sigTable[count,3] <- temp$K
  sigTable[count,4] <- temp$P
  count <- count + 1
}
```


Challenge 8: Anchovy phylogenetic signal
========================================================
type: prompt

Step 3: How'd we do? 

```r
reconAICcTable
```

```
      logbodysize trophic_position
BM       17.17923        -100.7962
OU       17.31522        -156.8098
EB       17.31522        -156.6843
trend    16.98079        -107.3646
white    49.82883        -159.0762
```

```r
sigTable
```

```
                       Lambda Lambda P-value         K K P-value
logbodysize      9.954688e-01   9.920151e-09 0.7099999     0.001
trophic_position 5.186526e-05   1.000000e+00 0.1009663     0.809
```

The Plan
====================================
type: section

- Introduction    
- Loading and Plotting Phylogenies
- Character Data
- Modeling Character Evolution
- Calculating Phylogenetic Signal
- **Ancestral Character State Reconstruction**


Ancestral Character State Reconstruction
====================================
type: section

Ancestral Character State Reconstruction
====================================

- The basic function
- How results are structured
- Continuous reconstruction
- Plot results
- Discrete reconstruction
- Plot results

Character reconstructions: Basic function
===============================================================


```r
recon <- ace(x = treeWData$data[,8], phy = treeWData$phy, 
            type = "continuous", method = "GLS", 
            corStruct = corBrownian(value=1, treeWData$phy))
```

Character reconstructions: How results are structured
===============================================================
- recon$ace numbers correspond to nodes in phylogeny
- recon$CI95 are confidence intervals for each node's reconstruction


```r
recon
```

```

    Ancestral Character Estimation

Call: ace(x = treeWData$data[, 8], phy = treeWData$phy, type = "continuous", 
    method = "GLS", corStruct = corBrownian(value = 1, treeWData$phy))

$ace
        12         13         14         15         16         17         18 
2589649968 1428149517 1272412103 1405633461 1393716796 1296733993 1235410112 
        19         20         21 
1084046522 1063892522  721486716 

$CI95
         [,1]       [,2]
12 2589649968 2589649968
13 1428149513 1428149522
14 1272412099 1272412107
15 1405633458 1405633465
16 1393716792 1393716799
17 1296733990 1296733997
18 1235410110 1235410114
19 1084046519 1084046525
20 1063892520 1063892525
21  721486714  721486719
```

Character reconstructions: Continuous reconstruction
===============================================================
The setup


```r
reconMatrix <- matrix(nrow=(length(treeWData$phy$tip.label) + treeWData$phy$Nnode), 
                      ncol=1)
colnames(reconMatrix) <- colnames(treeWData$data)[18]
rownames(reconMatrix) <- c(treeWData$phy$tip.label, 
                           seq(from = length(unlist(treeWData$phy$tip.label)) + 1, 
                               to = nrow(reconMatrix)))
```

Character reconstructions: Continuous reconstruction
===============================================================
Reconstruct niche breadth


```r
recon <- ace(treeWData$data[,8], treeWData$phy, type = "continuous", 
             method= "GLS", corStruct = corBrownian(value=1, treeWData$phy))
reconMatrix[,1] <- round(c(treeWData$data[,8], recon$ace), 2)
  #Putting data from reconstructions into a matrix.
reconMatrix <-as.matrix(reconMatrix)
```

Character reconstructions: Continuous reconstruction
===============================================================
Plotting reconstructed niche breadth


```r
library(RColorBrewer) # Accessory package with better color
colPal <- brewer.pal(9, "RdYlBu")
reconCol <- colPal[round((reconMatrix[,1]-min(reconMatrix[,1]))/(max(reconMatrix[,1])-min(reconMatrix[,1]))*8,0)+1]; #Determining what color will be displayed at each node and tip.
plot(treeWData$phy, main = "Reconstruction of niche breadth", label.offset = 2.5, no.margin = F) #Set down the plot.
nodelabels(text = rep(" ", 21-12+1), node = seq(12,21,1), cex = 0.65, frame = "circle", bg = reconCol[12:21]) #Label the nodes.
tiplabels(text = rep(" ", 11), seq(1,11,1), cex = 0.65, frame = "circle", bg = reconCol[1:11]) #Label the tips.

#Making a key for the plot
col.labels<-c("Narrow", "Broad")
color.legend(11,0,13,4,col.labels,colPal,gradient="y")
```

Character reconstructions: Continuous reconstruction
===============================================================
Plotting reconstructed niche breadth

<img src="IntroToCompPhy-figure/unnamed-chunk-51-1.png" title="plot of chunk unnamed-chunk-51" alt="plot of chunk unnamed-chunk-51" style="display: block; margin: auto;" />

Challenge 9: Reconstruct Anchovy body size
========================================================
type: prompt

Reconstruct and plot ancestral anchovy body size. Remember that the Ornstein Uhlenbeck evolutionary model was the best fit for this data. So the correlation structure you want to use is corMartins instead of corBrownian.

Challenge 9: Reconstruct Anchovy body size
========================================================
type: prompt

Step 1: Do the reconstruction.


```r
recon <- ace(chovyCharTab[,1], anchovyTreeWData$phy, 
             type = "continuous", 
             method= "GLS", corStruct = corBrownian(value=1, anchovyTreeWData$phy))
reconMatrix <- round(c(chovyCharTab[,1], recon$ace), 2)
  #Putting data from reconstructions into a matrix.
reconMatrix <-as.matrix(reconMatrix)
```

Challenge 9: Reconstruct Anchovy body size
========================================================
type: prompt

Step 2: Plot it.

```r
colPal <- brewer.pal(9, "PuOr")
reconCol <- colPal[round((reconMatrix[,1]-min(reconMatrix[,1]))/(max(reconMatrix[,1])-min(reconMatrix[,1]))*8,0)+1] #Determining what color will be displayed at each node and tip.

plot(anchovyTreeWData$phy, main = "Reconstruction of body size", label.offset = 2.5, no.margin = F) #Set down the plot.
nodelabels(text = rep(" ", 49), node = seq(51,99,1), cex = 0.65, frame = "circle", bg = reconCol[51:99]) #Label the nodes.
tiplabels(text = rep(" ", 50), seq(1,50,1), cex = 0.65, frame = "circle", bg = reconCol[1:50]) #Label the tips.

#Making a key for the plot
col.labels<-c("Small", "Large")
color.legend(0,50,10,40,col.labels,colPal,gradient="y")
```

Challenge 9: Reconstruct Anchovy body size
========================================================
type: prompt

Step 2: Plot it.  

<img src="IntroToCompPhy-figure/unnamed-chunk-54-1.png" title="plot of chunk unnamed-chunk-54" alt="plot of chunk unnamed-chunk-54" style="display: block; margin: auto;" />

Character reconstructions: Discrete reconstruction
===============================================================
The basic function


```r
discChar <- c(1, 1, 1, 3, 3, 1, 2, 2, 2, 2, 1)
  #Entering in discrete tip values
recon <- ace(discChar, treeWData$phy, type = "discrete", method = "ML", model = "ER");
  #Do reconstruction (equal likelihood for all transitions)
recon  #What it looks like
```

```

    Ancestral Character Estimation

Call: ace(x = discChar, phy = treeWData$phy, type = "discrete", method = "ML", 
    model = "ER")

    Log-likelihood: -7.502645 

Rate index matrix:
  1 2 3
1 . 1 1
2 1 . 1
3 1 1 .

Parameter estimates:
 rate index estimate std-err
          1    0.012  0.0081

Scaled likelihoods at the root (type '...$lik.anc' to get them for all nodes):
        1         2         3 
0.5242880 0.2443479 0.2313642 
```

Character reconstructions: Discrete reconstruction
===============================================================
What the ancestral state reconstructions look like.


```r
recon$lik.anc
```

```
                 1           2            3
 [1,] 0.5242879674 0.244347855 2.313642e-01
 [2,] 0.7988762539 0.129944244 7.117950e-02
 [3,] 0.9653742944 0.021355494 1.327021e-02
 [4,] 0.9927999747 0.004170209 3.029816e-03
 [5,] 0.6691495041 0.223404708 1.074458e-01
 [6,] 0.4347414214 0.386702471 1.785561e-01
 [7,] 0.0033931287 0.003106903 9.935000e-01
 [8,] 0.0292836926 0.957222775 1.349353e-02
 [9,] 0.0014809269 0.997760499 7.585745e-04
[10,] 0.0001420902 0.999761557 9.635278e-05
```

Character reconstructions: Discrete reconstruction
===============================================================
Plotting the reconstruction results.


```r
plot(treeWData$phy, main = "Reconstruction of preferred prey", cex.main = 5, label.offset = 2.5, no.margin = F, cex = 4)
  #Set down the plot
reconCol <- brewer.pal(3, "RdYlBu")
  #Setting up colors.
nodelabels(text = rep(" ", 21-12+1), node = seq(12,21,1), pie = recon$lik.anc, 
           piecol = c(reconCol[1], reconCol[2], reconCol[3]), cex = .9)
  #Label the nodes.
tiplabels(text = rep(" ", 11), seq(1,11,1), cex = 3, frame = "circle", bg = reconCol[discChar])
  #Label the tips.
legend(x = 0, y = 5, legend = c("Sardines", "Pompano", "Hamburgers"), fill = reconCol, bty = "n", cex = 5)
  #Plot the legend.
```

Character reconstructions: Discrete reconstruction
===============================================================
Plotting the reconstruction results.

<img src="IntroToCompPhy-figure/unnamed-chunk-58-1.png" title="plot of chunk unnamed-chunk-58" alt="plot of chunk unnamed-chunk-58" style="display: block; margin: auto;" />

Challenge 10: Reconstruct Anchovy diadromy
========================================================
type: prompt

Was the most recent common ancestor of anchovies more likely diadromous or not based on our data?

Hint: Use as.factor and as.numeric to transform your tip data from character strings to 1s and 2s!

Challenge 10: Reconstruct Anchovy diadromy
========================================================
type: prompt

Step 1: Do the reconstruction.


```r
#Set up the data
diadromy <- as.numeric(as.factor(anchovyTreeWData$data[,"Reg"]))
names(diadromy) <- rownames(anchovyTreeWData$data)

recon <- ace(diadromy, anchovyTreeWData$phy, type = "discrete", method = "ML", model = "ER")
  #Do reconstruction (equal likelihood for all transitions)
recon  #What it looks like
```

```

    Ancestral Character Estimation

Call: ace(x = diadromy, phy = anchovyTreeWData$phy, type = "discrete", 
    method = "ML", model = "ER")

    Log-likelihood: -23.89524 

Rate index matrix:
  1 2
1 . 1
2 1 .

Parameter estimates:
 rate index estimate std-err
          1   0.0056  0.0018

Scaled likelihoods at the root (type '...$lik.anc' to get them for all nodes):
        1         2 
0.0991219 0.9008781 
```

Challenge 10: Reconstruct Anchovy diadromy
========================================================
type: prompt

Step 2: Plot it.

Plotting the reconstruction results.


```r
plot(anchovyTreeWData$phy, main = "Evolution of Diadromy", cex.main = 5, label.offset = 2.5, no.margin = F, cex = 4)
  #Set down the plot
reconCol <- c("Red", "Blue")
  #Setting up colors.
nodelabels(text = rep(" ", 49), node = seq(51,99,1), 
           pie = recon$lik.anc, 
           piecol = c(reconCol[1], reconCol[2]), cex = .9)
  #Label the nodes.
tiplabels(text = rep(" ", 50), seq(1,50,1), 
          cex = 3, frame = "circle", bg = reconCol[diadromy])
  #Label the tips.
legend(x = -6, y = 50, legend = c("Diadromous", "Not"), 
       fill = reconCol, bty = "n", cex = 3)
  #Plot the legend.
```

Challenge 10: Reconstruct Anchovy diadromy
========================================================
type: prompt

Step 2: Plot it.

Plotting the reconstruction results.

<img src="IntroToCompPhy-figure/unnamed-chunk-61-1.png" title="plot of chunk unnamed-chunk-61" alt="plot of chunk unnamed-chunk-61" style="display: block; margin: auto;" />

Further Reading
====================================
type: section

Further Reading
========================================================

Analysis of Phylogenetics and Evolution with R
* By Emmanuel Paradis
  + Lead author of ape package
* http://ape-package.ird.fr/APER.html

Phytools
* Another CompPhy R package popularly used
* http://blog.phytools.org/
 
