
redNgray = function() c(red="#8B0000AA", gray="#D3D3D333")

# you can use mv.calout.detect on subsets of data defined by gene sets
# genesets = ivygapSE::makeGeneSets()

#' use shiny to explore univariate and multivariate outlier patterns in selected gene sets
#' @import shiny
#' @import parody
#' @importFrom beeswarm beeswarm
#' @importFrom MASS parcoord
#' @import maftools
#' @param genesets a named list of vectors of gene symbols
#' @examples
#' if (interactive()) conkApp()
#' @export
conkApp = function( genesets = conkout::glioSets47 ) {
mutatedGenes = c("default", sort(unique(conkout::tcga_gbm@data$Hugo_Symbol)))
#gopt = names(genesets)
#easy = substr(gopt,1,11)
#names(easy) = gopt
ui = fluidPage(
 sidebarLayout(
  sidebarPanel(
   helpText("conkout package outlier exploration, using TCGA-GBM hu133 expression data and glioblastoma-related gene sets from MSigDb"),
   selectInput("geneset", "geneset", names(genesets), names(genesets)[1]),
   uiOutput("picker"),
   numericInput("bp1", "bipl ax1", 1, min=1, max=10, step=1),
   numericInput("bp2", "bipl ax2", 2, min=2, max=10, step=1), 
#   selectInput("gene4mut", "gene for oncoprint", choices=mutatedGenes, selected="default"),
   width=4
   ),
  mainPanel(
   helpText("Tabs: beesOne for univariate, PCA for general princomp, biplot on PC selected as 'bipl ax', oncopl for maftools coOncoplot, which uses PoisonAlien TCGAMutations for GBM"),
   tabsetPanel(
    tabPanel("ParCo", plotOutput("curparco")),
    tabPanel("beesOne", plotOutput("beesOne")),
    tabPanel("PCA", plotOutput("prcomp")),
    tabPanel("biplot", plotOutput("bipl")),
    tabPanel("oncopl", plotOutput("oncop")),
    tabPanel("mutSumms", 
         textOutput("NumOLmsg"),
         plotOutput("oncopOL"), 
         textOutput("NumILmsg"),
         plotOutput("oncopIL"))
    )
  )
 )
)

server = function(input, output) {
  output$picker = renderUI(
     {
     curset = genesets[[input$geneset]]
     selectInput("curg", "gene for univ.", curset, curset[1])
     }
    )
  getCurrentOutliers = reactive({
     curgs = genesets[[input$geneset]]
     okg = intersect(rownames(conkout::gbmu133), curgs)
     texp = t(conkout::gbmu133[okg,])
     texp = texp[, order(colnames(texp))]
     list(olinds = mv.calout.detect(texp), texp=texp)
    })
  output$curparco = renderPlot({
     olstuff = getCurrentOutliers()
     colsToUse = rep(redNgray()[2], nrow(conkout::gbmu133))
     lwds = rep(1, nrow(conkout::gbmu133))
     if (sum(!is.na(curi <- olstuff$olinds$inds))>0)  {
        colsToUse[curi] = redNgray()[1]
        lwds[curi] = 3
        }
     par(las=2)
     parcoord(olstuff$texp, col=colsToUse, lwd=lwds,
       main=paste(input$geneset, " (Nmvout=", length(olstuff$olinds$inds), 
" Nin=", nrow(olstuff$texp)-length(olstuff$olinds$inds), ")", sep=""))
     })
   output$beesOne = renderPlot({
     olstuff = getCurrentOutliers()
     texp = olstuff$texp
     exps = texp[, input$curg]
     uniout = calout.detect(exps, method="GESD")
     mm = mean(exps)
     sd = sd(exps)
     dnlow = mm-2*sd
     dnhigh = mm+2*sd
     col2use = rep("gray", length(exps))
     col2usem = rep("gray", length(exps))
     col2useDN = rep("gray", length(exps))
     dnl = which(exps < dnlow)
     dnh = which(exps > dnhigh)
     if (length(dnl)>0) col2useDN[dnl] = "red"
     if (length(dnh)>0) col2useDN[dnh] = "red"
     if (!is.na(uniout$ind[1])) col2use[uniout$ind] = "red"
     if (!is.na(olstuff$olinds$ind[1])) col2usem[olstuff$olinds$ind] = "red"
     par(mfrow=c(1,3))
     beeswarm(exps, pwcol=col2useDN, pch=19, 
       main=paste0("DN 2SD rule: ", input$curg))
     beeswarm(exps, pwcol=col2use, pch=19, 
       main=paste0("univ GESD: ", input$curg))
     beeswarm(exps, pwcol=col2usem, pch=19, main=
      paste0("multivariate GESD:", input$curg))
     })
   output$prcomp = renderPlot({
     olstuff = getCurrentOutliers()
     texp = olstuff$texp
     pp = prcomp(texp)
     col2usem = rep(redNgray()[2], nrow(texp))
     if (!is.na(olstuff$olinds$ind[1])) col2usem[olstuff$olinds$ind] = redNgray()[1]
     pairs(pp$x[,1:4], pch=19, col=col2usem, cex=.9)
     })
   output$bipl = renderPlot({
     olstuff = getCurrentOutliers()
     texp = olstuff$texp
     pp = prcomp(texp)
     labs = rep(".", nrow(texp))
     if (!is.na(olstuff$olinds$ind[1])) labs[olstuff$olinds$ind] = "x"
     biplot(pp, xlabs=labs, choices=c(input$bp1, input$bp2))
     })
   prepMAF = reactive({
     olstuff = getCurrentOutliers()
     texp = olstuff$texp
     pp = prcomp(texp)
     labs = rep(".", nrow(texp))
     shn = substr(rownames(texp), 1, 12)
     if (!is.na(olstuff$olinds$ind[1])) {
       gr1 = shn[olstuff$olinds$ind] 
       mut1 = subsetMaf(conkout::tcga_gbm, gr1, mafObj=TRUE)
       gr2 = shn[-olstuff$olinds$ind] 
       mut2 = subsetMaf(conkout::tcga_gbm, gr2, mafObj=TRUE)
       }
     list(mutOL=mut1, mutIL=mut2) # outlier/inlier
   })
   output$oncop = renderPlot({
       mafstuff = prepMAF()
       coOncoplot(mafstuff$mutOL, mafstuff$mutIL) 
#       else if (input$gene4mut != "default") coOncoplot(mafstuff$mutOL, mafstuff$mutIL, genes=input$gene4mut) #, genes=input$geneset)
     })
   output$NumOLmsg = renderText({
       mafstuff = prepMAF()
       numOL = length(unique(mafstuff$mutOL@data$Tumor_Sample_Barcode))
       sprintf("MAF summary for %d outlying individuals with mutation info.", numOL)
       })
   output$NumILmsg = renderText({
       mafstuff = prepMAF()
       numIL = length(unique(mafstuff$mutIL@data$Tumor_Sample_Barcode))
       sprintf("MAF summary for %d inlying individuals with mutation info.", numIL)
       })
   output$oncopOL = renderPlot({
     mafstuff = prepMAF()
     plotmafSummary(mafstuff$mutOL)
     })
   output$oncopIL = renderPlot({
     mafstuff = prepMAF()
     plotmafSummary(mafstuff$mutIL)
     })
}
print(shinyApp(ui, server))
}
